#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "Rtypes.h"
#include "TFile.h"
#include "TTree.h"
#include "boost/expected/expected.hpp"

#include "analysis.h"
#include "branch_array_buffer_sizes.h"
#include "cutoffs.h"
#include "error.h"
#include "plotter.h"
#include "point.h"
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

namespace {

Expected<std::array<Int_t, 2>> FindStableLBBounds(const vector<string>* beam_mode) {
  auto this_func_name = "FindStableLBBounds";
  Int_t first_stable_LB = 0;
  for (const auto& status: *beam_mode) {
    if (status == "STABLE BEAMS") break;
    ++first_stable_LB;
  }

  Int_t last_stable_LB = beam_mode->size() - 1;
  for (auto status_it = beam_mode->rbegin();
       status_it != beam_mode->rend();
       ++status_it) {
    if (*status_it == "STABLE BEAMS") break;
    --last_stable_LB;
  }

  if (first_stable_LB == beam_mode->size()) {
    auto err_msg = "no stable beams found";
    return make_unexpected(make_shared<Error::Runtime>(err_msg, this_func_name));
  }
  else {
    assert(first_stable_LB <= last_stable_LB);
    return std::array<Int_t, 2>{first_stable_LB, last_stable_LB};
  }
}

}

SingleRunData::SingleRunData(std::string run_name, const Analysis& analysis)
  : run_name_(run_name),
    analysis_(analysis),
    nLB_(0),
    nCollisions_(0),
    LB_stability_offset_(0)
{
  if (!analysis_.use_start_of_fill_pedestals()){
    THROW_IF_ERR( ReadPedestals() )
  }
  THROW_IF_ERR( ReadTree() )
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
    unsigned iChannel = 0;
    for (const auto &channel: analysis_.channel_calibrations()) {
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

  auto LB_bounds = FindStableLBBounds(beam_mode);
  RETURN_IF_ERR( LB_bounds )
  auto first_stable_LB = LB_bounds->at(0);
  auto last_stable_LB = LB_bounds->at(1);
  // LB indexing starts at 1, I think
  LB_stability_offset_ = first_stable_LB + 1;

  if (analysis_.use_start_of_fill_pedestals()) {
    assert(pedestals_.size() == 0);
    auto iChannel = int{0};
    // Avoid reallocations by reusing the same vector
    vector<Float_t> pedestal_values;
    pedestal_values.reserve(first_stable_LB);
    for (const auto& channel: analysis_.channel_calibrations()) {
      auto channel_name = channel.first;
      for (Int_t iLB = 0; iLB < first_stable_LB; ++iLB) {
        if (lumi_BCM_temp[iLB] > gBCMLumiCutoff) {
          pedestal_values.push_back(currents_temp[iChannel][iLB]);
        }
      }
      auto pedestal_sum = std::accumulate(pedestal_values.begin(),
                                          pedestal_values.end(),
                                          0.0);
      auto pedestal_avg = pedestal_sum / pedestal_values.size();
      pedestals_.insert({channel_name, pedestal_avg});

      ++iChannel;
      pedestal_values.clear();
    }
  }

  if (analysis_.retrieve_currents()) {
    assert(currents_.size() == 0);

    auto iChannel = int{0};
    auto this_channel_currents = vector<Float_t>();
    this_channel_currents.reserve(last_stable_LB - first_stable_LB + 1);
    for (const auto& channel: pedestals_) {
      auto this_channel_name = channel.first;
      auto this_channel_pedestal = channel.second;

      for (int iLB = first_stable_LB; iLB <= last_stable_LB; ++iLB) {
        auto this_current = currents_temp[iChannel][iLB] -
                            this_channel_pedestal;
        this_channel_currents.push_back(this_current);
      }

      currents_.insert({std::move(this_channel_name),
                        std::move(this_channel_currents)});
      ++iChannel;
      this_channel_currents.clear();
    }
  }

  if (analysis_.retrieve_lumi_BCM()) {
    assert(lumi_BCM_.size() == 0);

    for (int iLB = first_stable_LB; iLB <= last_stable_LB; ++iLB) {
      lumi_BCM_.push_back( lumi_BCM_temp[iLB] );
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

// Calculates FCal luminosity from currents data for a run
Expected<Void> SingleRunData::CalcFCalLumi() {
  auto this_func_name = "SingleRunData::CalcFCalLumi";
  if (analysis_.verbose()) cout << "Calculating FCal luminosity" << endl;

  // Checks that these values were not previously calculated
  if (lumi_FCal_A_.size() > 0 || lumi_FCal_C_.size() > 0) {
    auto err_msg = "FCal lumi vector(s) already filled";
    return make_unexpected(make_shared<Error::Runtime>(err_msg, this_func_name));
  }

  // Checks that currents data exists
  if (currents_.size() == 0) {
    auto err_msg = "currents map empty";
    return make_unexpected(make_shared<Error::Runtime>(err_msg, this_func_name));
  }

  // All currents_ vectors should be the same size
  unsigned nLB_stable = currents_.begin()->second.size();
  if (nLB_stable == 0) {
    auto err_msg = "no stable lumi blocks";
    return make_unexpected(make_shared<Error::Runtime>(err_msg, this_func_name));
  }

  // Averages lumi measurement from all channels for each side and for each
  //   lumi block
  for (unsigned iLB = 0; iLB < nLB_stable; ++iLB) {
    Float_t lumi_A_temp = 0.0;
    Float_t lumi_C_temp = 0.0;
    unsigned num_channels_A = 0;
    unsigned num_channels_C = 0;

    // Calculates lumi for a channel using its current and calibration and
    //   adds to running total
    for (const auto& this_channel_currents: currents_) {
      const auto& channel_name = this_channel_currents.first;
      auto current = this_channel_currents.second.at(iLB);

      // Skip channel if current is ~0
      if (current < gFCalCurrentCutoff) continue;

      // Perf sore point? (like it matters)
      auto intercept = analysis_.channel_calibrations().at(channel_name).intercept;
      auto slope = analysis_.channel_calibrations().at(channel_name).slope;

      if (channel_name.at(1) == '1') { // channel from C-side
        lumi_C_temp += current*slope + intercept;
        ++num_channels_C;
      }
      else {
        lumi_A_temp += current*slope + intercept;
        ++num_channels_A;
      }
    }

    // Computes the average lumi for this lumi block and adds to lumi vector.
    // Adds 0 if no channels had high enough current
    if (num_channels_A == 0) {
      lumi_FCal_A_.push_back(0.0);
    }
    else {
      lumi_FCal_A_.push_back(analysis_.corr_A()*lumi_A_temp / num_channels_A);
    }
    if (num_channels_C == 0) {
      lumi_FCal_C_.push_back(0.0);
    }
    else {
      lumi_FCal_C_.push_back(analysis_.corr_C()*lumi_C_temp / num_channels_C);
    }
  }

  return Void();
}

// Calculates <mu> from FCal luminosity for a run
Expected<Void> SingleRunData::CalcFCalMu() {
  auto this_func_name = "SingleRunData::CalcFCalMu";

  if (analysis_.verbose()) cout << "Calculating FCal <mu>" << endl;

  // Checks that FCal lumi data exists
  if (lumi_FCal_A_.size() == 0 && lumi_FCal_C_.size() == 0) {
    auto err_msg = "FCal lumi vector(s) have not been filled";
    return make_unexpected(make_shared<Error::Runtime>(err_msg, this_func_name));
  }
  // Checks that mu data has not already been calculated
  if (mu_FCal_A_.size() > 0 || mu_FCal_C_.size() > 0) {
    auto err_msg = "FCal mu vector(s) already filled";
    return make_unexpected(make_shared<Error::Runtime>(err_msg, this_func_name));
  }
  // Checks that value for number of collisions is nonzero
  if (nCollisions_ == 0) {
    auto err_msg = "num_collisions == 0";
    return make_unexpected(make_shared<Error::Runtime>(err_msg, this_func_name));
  }

  Float_t conversion_factor = analysis_.x_sec() / (nCollisions_ * analysis_.f_rev());
  // Cross-section conversion factor from 7 -> 8 TeV
  if (run_name_.at(0) == '2') conversion_factor /= 1.05;

  // Calculates <mu> for each lumi block
  for (const auto &lumi: lumi_FCal_A_) {
    mu_FCal_A_.push_back(lumi*conversion_factor);
  }
  for (const auto &lumi: lumi_FCal_C_) {
    mu_FCal_C_.push_back(lumi*conversion_factor);
  }
  return Void();
}

Expected<Void> SingleRunData::CreateLumiCurrentPlots() const {
  const char* this_func_name = "SingleRunData::CreateLumiCurrentPlots";
  if (analysis_.verbose()) cout << "    " << "Making lumi vs. current plots" << endl;

  LumiCurrentPlotOptions plot_options(analysis_.params_filepath());
  plot_options.set_run_name(run_name_);
  std::map<string, FitResults> fit_results;

  if (plot_options.do_individual()) {
    for (const auto &channel: currents_) {
      auto this_channel_name = channel.first;
      plot_options.set_channel_name(this_channel_name);

      auto lumi_current_points = PointVectorFromVectors(
          lumi_BCM_,
          channel.second);

      CONTINUE_IF_ERR(lumi_current_points);

      auto this_channel_fit_results =
        Plotter::PlotLumiCurrent(*lumi_current_points,
                                 plot_options);

      CONTINUE_IF_ERR(this_channel_fit_results)

      if (pedestals_.at(this_channel_name) > 20.0) {
        this_channel_fit_results->is_short = true;
      }
      else {
        this_channel_fit_results->is_short = false;
      }

      if (plot_options.do_fit()) {
        fit_results.insert({this_channel_name, *this_channel_fit_results});
      }
    } //channels loop
  }

  if (plot_options.do_sum()) {
    vector<Float_t> channel_currents_sum_A;
    vector<Float_t> channel_currents_sum_C;
    for (const auto &channel: analysis_.channel_calibrations()) {
      auto this_channel_name = channel.first;
      auto this_channel_current = currents_.at(this_channel_name);

      vector<Float_t> *this_side_channel_currents_sum;
      if (this_channel_name.at(1) == '1') {
        this_side_channel_currents_sum = &channel_currents_sum_A;
      }
      else if (this_channel_name.at(1) == '8') {
        this_side_channel_currents_sum = &channel_currents_sum_C;
      }
      else {
        return make_unexpected(make_shared<Error::Logic>("Invalid channel name", this_func_name));
      }

      if (this_side_channel_currents_sum->size() == 0) {
        *this_side_channel_currents_sum = this_channel_current;
      }
      else {
        std::transform(this_side_channel_currents_sum->begin(),
                       this_side_channel_currents_sum->end(),
                       this_channel_current.begin(),
                       this_side_channel_currents_sum->begin(),
                       std::plus<Float_t>());
      }
    }

    // TODO: figure out how to deal with sums
    auto points_sum_A = PointVectorFromVectors(lumi_BCM_,
                                               channel_currents_sum_A);
    plot_options.set_channel_name("Sum_A");
    if (points_sum_A.valid()) {
      TRY( Plotter::PlotLumiCurrent(*points_sum_A, plot_options) )
    }

    auto points_sum_C = PointVectorFromVectors(lumi_BCM_,
                                               channel_currents_sum_C);
    plot_options.set_channel_name("Sum_C");
    if (points_sum_C.valid()) {
      TRY( Plotter::PlotLumiCurrent(*points_sum_C, plot_options) )
    }
  }

  if (plot_options.do_fit()) {
    LOG_IF_ERR( Plotter::WriteFitResultsToTree(fit_results,
                                   plot_options) )

    LOG_IF_ERR( Plotter::WriteCalibrationToText(fit_results,
                                    plot_options) );

    LOG_IF_ERR( Plotter::GeometricAnalysisOfFitResults(fit_results,
                                           plot_options) );
  }

  return Void();
}

// Create those plots which use data from only a single sample
Expected<Void> SingleRunData::CreateSingleRunPlots() const {
  for (const auto &plot_type: analysis_.plot_types()) {
    if (plot_type == "lumi_current") {
      if (analysis_.verbose()) {
        cout << "Creating plots for sample " << run_name_ << endl;
      }
      TRY( CreateLumiCurrentPlots() )
    }
  }
  // Remove the TRY above if other single run plots are added
  return Void();
}

