#include <algorithm>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "Rtypes.h"
#include "TFile.h"
#include "TTree.h"
#include "boost/expected/expected.hpp"

#include "analysis.h"
#include "branch_array_buffer_sizes.h"
#include "error.h"
#include "util.h"
#include "single_run_data.h"

using std::string;
using std::vector;

using std::cout;
using std::cerr;
using std::endl;

using std::make_shared;

using boost::make_unexpected;

using Error::Expected;

SingleRunData::SingleRunData(std::string run_name, const Analysis& analysis)
  : run_name_(run_name),
    analysis_(analysis),
    nLB_(0),
    nCollisions_(0),
    LB_stability_offset_(0)
{
  THROW_IF_ERR( ReadPedestals() )
  THROW_IF_ERR( ReadTree() )
}

void SingleRunData::InitCurrentsMap()
{
  for (const auto& run: pedestals_) {
    const auto& channel_name = run.first;
    currents_.insert(std::make_pair(channel_name, vector<Float_t>()));
  }
}

// Reads in pedestal values for each of the channels specified in analysis_
Expected<Void> SingleRunData::ReadPedestals()
{
  auto this_func_name = "SingleRunData::ReadPedestals";

  auto pedestals_filepath = analysis_.pedestals_dir() + run_name_ + ".dat";
  std::ifstream pedestals_file(pedestals_filepath);
  if (!pedestals_file) {
    return make_unexpected(make_shared<Error::File>(pedestals_filepath, this_func_name));
  }

  string channel_name;
  Float_t pedestal;
  Float_t rms_noise;
  while (pedestals_file >> channel_name >> pedestal >> rms_noise) {
    const auto& channels_list = analysis_.channel_calibrations();
    if ( std::find_if(channels_list.begin(), channels_list.end(),
                      [&channel_name](const auto& run) {
                        return channel_name == run.first;
                      })
        != channels_list.end() ) {
      pedestals_.insert({channel_name, pedestal});
    }
  }

  if (analysis_.channel_calibrations().size() != pedestals_.size()) {
    auto err_msg = "missing pedestal for a channel";
    return make_unexpected(make_shared<Error::Runtime>(err_msg, this_func_name));
  }

  return Void();
}

// Reads relevant data in from the file {analysis_.trees_dir()}/{run_name_}.root
//   and copies it to member variables.
Expected<Void> SingleRunData::ReadTree()
{
  const char *this_func_name = "SingleRunData::ReadTree";

  if (analysis_.verbose()) cout << "Analysing sample " << run_name_ << endl;

  string filepath = analysis_.trees_dir()+run_name_+".root";
  TFile *this_file = TFile::Open(filepath.c_str());
  if (!this_file) {
    return make_unexpected(make_shared<Error::File>(filepath, this_func_name));
  }

  string tree_name = "t_"+run_name_;
  TTree *this_tree = static_cast<TTree*>( this_file->Get(tree_name.c_str()) );
  if (!this_tree) {
    auto err_msg = "tree `"+tree_name+"` not found in file `"+filepath+"`";
    return make_unexpected(make_shared<Error::Runtime>(err_msg, this_func_name));
  }

  this_tree->SetBranchAddress("StartOfRun", &timestamp_);
  this_tree->SetBranchAddress("NLB", &nLB_);
  this_tree->SetBranchAddress("ncoll", &nCollisions_);

  // Allocate buffers big enough to hold the data (can't use dynamic memory
  //   because of how TTrees work)
  Float_t currents_temp[gMaxNumChannels][gMaxNumLB];
  Float_t lumi_BCM_temp[gMaxNumLB];
  if (analysis_.retrieve_currents()) {
    InitCurrentsMap();
    unsigned iChannel = 0;
    for (const auto &channel: currents_) {
      auto branch_name = "current_"+channel.first;
      this_tree->SetBranchAddress((branch_name).c_str(), &currents_temp[iChannel]);
      ++iChannel;
    }
  }
  if (analysis_.retrieve_lumi_BCM()) {
    this_tree->SetBranchAddress("ofl_lumi_pref",&lumi_BCM_temp);
  }

  int quality[gMaxNumLB];
  vector<string> *beam_mode = nullptr;
  this_tree->SetBranchAddress("quality", &quality);
  this_tree->SetBranchAddress("mode", &beam_mode);

  // Populate those variables which have been set to branches
  this_tree->GetEntry(0);
  this_file->Close();

  HardcodenLBIfMissingFromTree();

  for (int iLB = 0; iLB < nLB_; ++iLB) {
    if (quality[iLB] && beam_mode->at(iLB) == "STABLE BEAMS") {
      if (LB_stability_offset_ == 0) {
        LB_stability_offset_ = iLB + 1; //LB numbers are 1-indexed (I think)
      }
      // Branch prediction hopefully removes the inefficiency here
      if (analysis_.retrieve_currents()) {
        unsigned iChannel = 0;
        for (const auto& channel: pedestals_) {
          string this_channel_name = channel.first;
          Float_t this_channel_pedestal = channel.second;

          auto this_current = currents_temp[iChannel][iLB] -
                              this_channel_pedestal;
          currents_.at(this_channel_name).push_back(this_current);
          ++iChannel;
        }
      }
      if (analysis_.retrieve_lumi_BCM()) {
        lumi_BCM_.push_back( lumi_BCM_temp[iLB] );
      }
    }
  }

  return Void();
}

Expected<Void> SingleRunData::CreateBenedettoOutput() const
{
  auto this_func_name = "SingleRunData::CreateBenedettoOutput";

  auto output_dir = analysis_.benedetto_output_dir();
  TRY( Util::mkdir(output_dir) )
  auto output_filepath = output_dir+run_name_+".dat";
  std::ofstream output_file(output_filepath);
  if (!output_file.is_open()) {
    return make_unexpected(make_shared<Error::File>(output_filepath, this_func_name));
  }

  auto num_points = mu_FCal_A_.size();
  for (unsigned iLB = 0; iLB < num_points; ++iLB) {
    output_file << (iLB + LB_stability_offset_) << ' ' << mu_FCal_A_.at(iLB) << ' '
                << mu_FCal_C_.at(iLB) << std::endl;
  }
  return Void();
}

const std::vector<std::pair<string, Int_t>> MISSING_nLB { {"206955", 1368},
                                                          {"208642", 465},
                                                          {"211620", 801} };
void SingleRunData::HardcodenLBIfMissingFromTree()
{
  auto run_ptr = std::find_if(MISSING_nLB.begin(), MISSING_nLB.end(),
                              [this](const auto& run) {
                                return run_name_ == run.first;
                              });
  if (run_ptr != MISSING_nLB.end()) {
    if (analysis_.verbose()) {
      cout << "Manually setting nColl for run " << run_name_ << endl;
    }
    nLB_ = run_ptr->second;
  }
}
