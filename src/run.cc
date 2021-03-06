#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "Rtypes.h"
#include "TFile.h"
#include "TTree.h"

#include "boost/expected/expected.hpp"
#include "boost/container/flat_map.hpp"

#include "analysis.h"
#include "branch_array_buffer_sizes.h"
#include "cutoffs.h"
#include "error.h"
#include "fcal_region.h"
#include "plotter.h"
#include "point.h"
#include "util.h"
#include "run.h"

using std::string;
using std::vector;

using std::cout;
using std::cerr;
using std::endl;

template<typename K, typename V>
using map = boost::container::flat_map<K, V>;

using std::make_shared;

using boost::make_unexpected;

using Error::Expected;

namespace {

void DumpLArAndOflLumiComparison(const Run& run)
{
  const auto& lumi_ofl = run.lumi_ofl();
  const auto& lumi_LAr_A = run.lumi_LAr_A();
  const auto& lumi_LAr_C = run.lumi_LAr_C();

  cout << lumi_ofl.size() << endl;
  cout << lumi_LAr_A.size() << endl;
  cout << lumi_LAr_C.size() << endl;

  assert(lumi_ofl.size() == lumi_LAr_A.size());
  assert(lumi_ofl.size() == lumi_LAr_C.size());

  const auto n_LB = lumi_ofl.size();
  double A_ratio_avg = 0.;
  double C_ratio_avg = 0.;
  int n_LB_used_in_A_avg = 0;
  int n_LB_used_in_C_avg = 0;

  for (int i = 0; i < n_LB; ++i) {
    auto this_LB_lumi_ofl = lumi_ofl.at(i);
    auto A_diff = std::abs(lumi_LAr_A.at(i) - this_LB_lumi_ofl);
    auto A_ratio = A_diff / this_LB_lumi_ofl * 100;
    auto C_diff = std::abs(lumi_LAr_C.at(i) - this_LB_lumi_ofl);
    auto C_ratio = C_diff / this_LB_lumi_ofl * 100;
    auto width = 8;
    cout << std::fixed << std::setw(4) << std::setprecision(0) << i+run.LB_stability_offset() << ' '
                       << std::setw(width) << std::setprecision(2) << this_LB_lumi_ofl << ' '
                       << std::setw(width) << std::setprecision(4) << A_ratio << ' '
                       << std::setw(width) << std::setprecision(4) << C_ratio << ' ' << endl;
    if (A_ratio < 5) {
      A_ratio_avg += A_ratio;
      ++n_LB_used_in_A_avg;
    }
    if (C_ratio < 5) {
      C_ratio_avg += C_ratio;
      ++n_LB_used_in_C_avg;
    }
  }
  cout << "A ratio average (excluding extreme outliers): "
       << A_ratio_avg/n_LB_used_in_A_avg << endl;
  cout << "C ratio average (excluding extreme outliers): "
       << C_ratio_avg/n_LB_used_in_C_avg << endl;
}

Expected<std::array<Int_t, 2>> FindStableLBBounds(const vector<string>* beam_mode)
{
  auto this_func_name = "FindStableLBBounds";

  if (beam_mode->size() == 0) {
    auto err_msg = "empty beam_mode vector";
    return make_unexpected(make_shared<Error::Runtime>(err_msg, this_func_name));
  }

  Int_t first_stable_LB = 0;
  for (const auto& status: *beam_mode) {
    if (status == "STABLE BEAMS") break;
    ++first_stable_LB;
  }

  if (first_stable_LB == beam_mode->size()) {
    auto err_msg = "no stable beams found";
    return make_unexpected(make_shared<Error::Runtime>(err_msg, this_func_name));
  }

  Int_t last_stable_LB = beam_mode->size() - 1;
  for (auto status_it = beam_mode->rbegin();
       status_it != beam_mode->rend();
       ++status_it) {
    if (*status_it == "STABLE BEAMS") break;
    --last_stable_LB;
  }

  assert(first_stable_LB <= last_stable_LB);
  return std::array<Int_t, 2>{first_stable_LB, last_stable_LB};
}

Expected<std::array<Int_t, 2>> FindReadyForPhysicsLBBounds(const vector<int>& quality)
{
  auto this_func_name = "FindReadyForPhysicsLBBounds";

  if (quality.size() == 0) {
    auto err_msg = "empty quality vector";
    return make_unexpected(make_shared<Error::Runtime>(err_msg, this_func_name));
  }

  Int_t first_stable_LB = 0;
  for (const auto& status: quality) {
    if (status == 1) break;
    ++first_stable_LB;
  }

  if (first_stable_LB == quality.size()) {
    auto err_msg = "no Ready-for-Physics LBs found";
    return make_unexpected(make_shared<Error::Runtime>(err_msg, this_func_name));
  }

  Int_t last_stable_LB = quality.size() - 1;
  for (auto status_it = quality.rbegin();
       status_it != quality.rend();
       ++status_it) {
    if (*status_it == 1) break;
    --last_stable_LB;
  }

  assert(first_stable_LB <= last_stable_LB);
  return std::array<Int_t, 2>{first_stable_LB, last_stable_LB};
}

template<class T>
vector<T> CArrayToVec(T* c_array, size_t size)
{
  vector<T> vec;
  vec.reserve(size);
  for (auto i = 0; i < size; ++i) {
    vec.push_back(c_array[i]);
  }
  return vec;
}

Expected<Int_t> FindStartOfAdjustLB(const vector<string>* beam_mode)
{
  auto this_func_name = "FindStartOfAdjustLB";

  if (beam_mode->size() == 0) {
    auto err_msg = "empty beam_mode vector";
    return make_unexpected(make_shared<Error::Runtime>(err_msg, this_func_name));
  }

  auto initial_beam_mode = beam_mode->at(0);
  if (initial_beam_mode == "STABLE BEAMS" ||
      initial_beam_mode == "ADJUST") {
    auto err_msg = "stable or adjust beams at first lumi block";
    return make_unexpected(make_shared<Error::Runtime>(err_msg, this_func_name));
  }

  Int_t start_of_adjust_LB = 0;
  for (const auto& status: *beam_mode) {
    if (status == "ADJUST") break;
    ++start_of_adjust_LB;
  }

  return start_of_adjust_LB;
}

float CalculateAvgBeamspotZ(const vector<float>& beamspot_z)
{
  float sum = 0;
  int n = 0;
  for (auto this_beamspot_z: beamspot_z) {
    if (this_beamspot_z == Run::BeamspotPlaceholderVal()) {
      continue;
    }
    else {
      sum += this_beamspot_z;
      ++n;
    }
  }

  if (n == 0) {
    return Run::BeamspotPlaceholderVal();
  }
  else {
    return sum/n;
  }
}

}

Run::Run(std::string run_name, Analysis* analysis)
  : run_name_(std::move(run_name)),
    analysis_(analysis)
{}

Expected<Void> Run::Init()
{
  if (!analysis_->use_start_of_fill_pedestals()){
    TRY( ReadPedestals() )
  }
  if (analysis_->use_baseline_subtraction_from_fit()) {
    channel_calibrations_ = analysis_->channel_calibrations();
    auto calibrations_filepath = analysis_->calibrations_dir()+run_name_+".dat";
    TRY( Analysis::ReadCalibrations(&channel_calibrations_, calibrations_filepath) )
  }
  TRY( ReadTree() )
  return Void();
}

// Reads in pedestal values for each of the channels specified in analysis_
Expected<Void> Run::ReadPedestals()
{
  auto this_func_name = "Run::ReadPedestals";

  auto pedestals_filepath = analysis_->pedestals_dir() + run_name_ + ".dat";
  std::ifstream pedestals_file(pedestals_filepath);
  if (!pedestals_file) {
    return make_unexpected(make_shared<Error::File>(pedestals_filepath, this_func_name));
  }

  string channel_name;
  Float_t pedestal;
  Float_t rms_noise;
  while (pedestals_file >> channel_name >> pedestal >> rms_noise) {
    const auto& channels_list = analysis_->channel_calibrations();
    if ( std::find_if(channels_list.begin(), channels_list.end(),
                      [&channel_name](const auto& run) {
                        return channel_name == run.first;
                      })
        != channels_list.end() ) {
      pedestals_.insert({{channel_name, pedestal}});
    }
  }

  if (analysis_->channel_calibrations().size() != pedestals_.size()) {
    auto err_msg = "missing pedestal for channel " + run_name_;
    return make_unexpected(make_shared<Error::Runtime>(err_msg, this_func_name));
  }

  return Void();
}

// Reads relevant data in from the file {analysis_->trees_dir()}/{run_name_}.root
//   and copies it to member variables.
Expected<Void> Run::ReadTree()
{
  auto this_func_name = "Run::ReadTree";

  if (analysis_->verbose()) cout << "Analysing sample " << run_name_ << endl;

  string luminosity_file_path = analysis_->trees_dir()+run_name_+".root";
  TFile *luminosity_file = TFile::Open(luminosity_file_path.c_str());
  if (!luminosity_file) {
    return make_unexpected(make_shared<Error::File>(luminosity_file_path, this_func_name));
  }

  string tree_name = "t_"+run_name_;
  TTree *luminosity_tree = static_cast<TTree*>( luminosity_file->Get(tree_name.c_str()) );
  if (!luminosity_tree) {
    auto err_msg = "tree `"+tree_name+"` not found in file `"+luminosity_file_path+"`";
    return make_unexpected(make_shared<Error::Runtime>(err_msg, this_func_name));
  }

  luminosity_tree->SetBranchAddress("start_of_run", &timestamp_);
  luminosity_tree->SetBranchAddress("n_LB", &n_LB_);

  Int_t n_bunches_tmp;
  luminosity_tree->SetBranchAddress("n_bunches", &n_bunches_tmp);

  // Allocate buffers big enough to hold the data (can't use dynamic memory
  //   because of how TTrees work (I think))
  Float_t currents_temp[gMaxNumChannels][gMaxNumLB];
  Float_t lumi_ofl_temp[gMaxNumLB];
  if (analysis_->retrieve_currents()) {
    string current_file_path = analysis_->currents_dir()+run_name_+".root";
    TFile *current_file = TFile::Open(current_file_path.c_str());
    if (!current_file) {
      return make_unexpected(make_shared<Error::File>(current_file_path, this_func_name));
    }

    string tree_name = "t_run_"+run_name_+"_currents";
    TTree *current_tree = static_cast<TTree*>( current_file->Get(tree_name.c_str()) );
    if (!current_tree) {
      auto err_msg = "tree `"+tree_name+"` not found in file `"+current_file_path+"`";
      return make_unexpected(make_shared<Error::Runtime>(err_msg, this_func_name));
    }
    for (auto i_LB = 0; i_LB < current_tree->GetEntries(); ++i_LB) {
      unsigned iChannel = 0;
      for (const auto &channel: analysis_->channel_calibrations()) {
        auto branch_name = channel.first;
        current_tree->SetBranchAddress((branch_name).c_str(), &currents_temp[iChannel][i_LB]);
        ++iChannel;
      }
      current_tree->GetEntry(i_LB);
    }
    current_file->Close();
  }
  if (analysis_->retrieve_lumi_ofl()) {
    luminosity_tree->SetBranchAddress(analysis_->reference_lumi_algo().c_str(), &lumi_ofl_temp);
  }

  // Ready for physics flag state for each LB
  Int_t RFP_flag_arr[gMaxNumLB];
  Float_t beamspot_z_arr[gMaxNumLB];
  // Ramp, Stable, etc.
  vector<string> *beam_mode = nullptr;
  Float_t beamspot_z[gMaxNumLB];
  luminosity_tree->SetBranchAddress("LB_quality", &RFP_flag_arr);
  luminosity_tree->SetBranchAddress("beam_mode", &beam_mode);

  if (analysis_->retrieve_beamspot()) {
    luminosity_tree->SetBranchAddress("beamspot_z", &beamspot_z);
    luminosity_tree->SetBranchAddress("beamspot_z_avg", &avg_beamspot_z_);
    beamspot_z_ = CArrayToVec<Float_t>(beamspot_z_arr, n_LB_);
  }

  // Populate those variables which have been set to branches
  luminosity_tree->GetEntry(0);
  luminosity_file->Close();

  n_bunches_ = n_bunches_tmp;

  if (analysis_->retrieve_mu_LAr() && !analysis_->n_bunches_from_file().empty()) GetExternalNBunches();

  RFP_flag_ = CArrayToVec<Int_t>(RFP_flag_arr, n_LB_);

  auto LB_bounds = GetLBBounds();
  RETURN_IF_ERR( LB_bounds )
  auto lower_LB_bound = LB_bounds->at(0);
  auto upper_LB_bound = LB_bounds->at(1);
  // LB indexing starts at 1 within ATLAS
  LB_stability_offset_ = lower_LB_bound + 1;

  if (analysis_->use_start_of_fill_pedestals()) {
    assert(pedestals_.size() == 0);

    auto start_of_adjust_LB = FindStartOfAdjustLB(beam_mode);

    RETURN_IF_ERR( start_of_adjust_LB )

    auto iChannel = 0;
    // Avoid reallocations by reusing the same vector
    // I could do this sans vector, but idc about the memory usage of a < 100-
    //   element vector
    vector<Float_t> pedestal_values;
    pedestal_values.reserve(*start_of_adjust_LB);
    for (const auto& channel: analysis_->channel_calibrations()) {
      auto channel_name = channel.first;
      for (Int_t i_LB = 0; i_LB < lower_LB_bound; ++i_LB) {
        //if (channel_name == "M80C0") cout << lumi_ofl_temp[i_LB] << endl;
        // We want LBs from the start of the fill before beam align
        auto this_channel_LB_current = currents_temp[iChannel][i_LB];

        // Sometimes we get current of exactly 0 at the start of beam inject or
        // setup. ofl offline lumi should always be ~0 before beam adjust
        if (lumi_ofl_temp[i_LB] < gOflLumiCutoff &&
            this_channel_LB_current > gEpsilon) {
          pedestal_values.push_back(this_channel_LB_current);
        }
      }
      // Arbitrary restriction on how many good pedestal readings there must be
      const auto min_pedestal_values = 3;
      if (pedestal_values.size() < min_pedestal_values) {
        auto err_msg = "insufficient currents (< " +
                       std::to_string(min_pedestal_values) +
                       ") in channel " + channel_name +
                       " to form the avg pedestal calculation";
        return make_unexpected(make_shared<Error::Runtime>(err_msg, this_func_name));
      }

      auto pedestal_sum = std::accumulate(pedestal_values.begin(),
                                          pedestal_values.end(),
                                          0.0);
      auto pedestal_avg = pedestal_sum / pedestal_values.size();
      pedestals_.insert({{channel_name, pedestal_avg}});

      ++iChannel;
      pedestal_values.clear();
    }
  }

  if (analysis_->retrieve_currents()) {
    assert(currents_.size() == 0);

    auto iChannel = 0;
    auto this_channel_currents = vector<Float_t>();
    this_channel_currents.reserve(upper_LB_bound - lower_LB_bound + 1);
    for (const auto& channel: pedestals_) {
      auto this_channel_name = channel.first;
      auto this_channel_pedestal = channel.second;
      if (analysis_->use_baseline_subtraction_from_fit()) {
        this_channel_pedestal -= channel_calibrations_.at(this_channel_name).intercept/channel_calibrations_.at(this_channel_name).slope;
      }

      for (int i_LB = lower_LB_bound; i_LB <= upper_LB_bound; ++i_LB) {
        auto this_current = currents_temp[iChannel][i_LB] -
                            this_channel_pedestal;
        this_channel_currents.push_back(this_current);
      }

      currents_.insert({{std::move(this_channel_name),
                         std::move(this_channel_currents)}});
      ++iChannel;
      this_channel_currents.clear();
    }
  }

  if (analysis_->retrieve_lumi_ofl()) {
    assert(lumi_ofl_.size() == 0);

    // Mu-dependence correction evaluated by Vincent Hedberg and shown in this presentation by
    // Benedetto Giacobbe:
    // https://indico.cern.ch/event/435624/contributions/2209977/attachments/1293235/1927107/LTF_16062016.pdf
    if (analysis_->apply_LUCID_mu_corr()) {
      auto conversion_factor = analysis_->x_sec() / (n_bunches() * analysis_->f_rev());
      for (int i_LB = lower_LB_bound; i_LB <= upper_LB_bound; ++i_LB) {
        auto lumi = lumi_ofl_temp[i_LB];
        // have to convert to mu, apply the correction, then convert back to lumi
        auto mu = lumi*conversion_factor;
        auto mu_corr = -0.002*mu*mu + 1.008*mu; // this is the correction
        auto lumi_corr = mu_corr/conversion_factor;
        lumi_ofl_.push_back( lumi_corr );
      }
    }
    else {
      for (int i_LB = lower_LB_bound; i_LB <= upper_LB_bound; ++i_LB) {
        lumi_ofl_.push_back( lumi_ofl_temp[i_LB] );
      }
    }
  }

  return Void();
}

// Reads relevant data in from the file {analysis_->trees_dir()}/{run_name_}.root
//   and copies it to member variables.
Expected<Void> Run::ReadTreeLegacy()
{
  auto this_func_name = "Run::ReadTree";

  if (analysis_->verbose()) cout << "Analysing sample " << run_name_ << endl;

  string filepath = analysis_->trees_dir()+run_name_+".root";
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
  this_tree->SetBranchAddress("NLB", &n_LB_);
  this_tree->SetBranchAddress("ncoll", &n_bunches_);

  // Allocate buffers big enough to hold the data (can't use dynamic memory
  //   because of how TTrees work (I think))
  Float_t currents_temp[gMaxNumChannels][gMaxNumLB];
  Float_t lumi_ofl_temp[gMaxNumLB];
  if (analysis_->retrieve_currents()) {
    unsigned iChannel = 0;
    for (const auto &channel: analysis_->channel_calibrations()) {
      auto branch_name = "current_"+channel.first;
      this_tree->SetBranchAddress((branch_name).c_str(), &currents_temp[iChannel]);
      ++iChannel;
    }
  }
  if (analysis_->retrieve_lumi_ofl()) {
    this_tree->SetBranchAddress("ofl_lumi_pref",&lumi_ofl_temp);
  }

  // Ready for physics flag state for each LB
  Int_t RFP_flag_arr[gMaxNumLB];
  Float_t beamspot_z_arr[gMaxNumLB];
  // Ramp, Stable, etc.
  vector<string> *beam_mode = nullptr;
  Float_t beamspot_z[gMaxNumLB];
  this_tree->SetBranchAddress("quality", &RFP_flag_arr);
  this_tree->SetBranchAddress("mode", &beam_mode);

  if (analysis_->retrieve_beamspot()) {
    this_tree->SetBranchAddress("beamspot_z", &beamspot_z);
    this_tree->SetBranchAddress("beamspot_z_avg", &avg_beamspot_z_);
    beamspot_z_ = CArrayToVec<Float_t>(beamspot_z_arr, n_LB_);
  }

  // Populate those variables which have been set to branches
  this_tree->GetEntry(0);
  this_file->Close();

  GetExternalNBunches();

  RFP_flag_ = CArrayToVec<Int_t>(RFP_flag_arr, n_LB_);

  auto LB_bounds = GetLBBounds();
  RETURN_IF_ERR( LB_bounds )
  auto lower_LB_bound = LB_bounds->at(0);
  auto upper_LB_bound = LB_bounds->at(1);
  // LB indexing starts at 1 within ATLAS
  LB_stability_offset_ = lower_LB_bound + 1;

  if (analysis_->use_start_of_fill_pedestals()) {
    assert(pedestals_.size() == 0);

    auto start_of_adjust_LB = FindStartOfAdjustLB(beam_mode);

    RETURN_IF_ERR( start_of_adjust_LB )

    auto iChannel = 0;
    // Avoid reallocations by reusing the same vector
    // I could do this sans vector, but idc about the memory usage of a < 100-
    //   element vector
    vector<Float_t> pedestal_values;
    pedestal_values.reserve(*start_of_adjust_LB);
    for (const auto& channel: analysis_->channel_calibrations()) {
      auto channel_name = channel.first;
      for (Int_t i_LB = 0; i_LB < lower_LB_bound; ++i_LB) {
        //if (channel_name == "M80C0") cout << lumi_ofl_temp[i_LB] << endl;
        // We want LBs from the start of the fill before beam align
        auto this_channel_LB_current = currents_temp[iChannel][i_LB];

        // Sometimes we get current of exactly 0 at the start of beam inject or
        // setup. ofl offline lumi should always be ~0 before beam adjust
        if (lumi_ofl_temp[i_LB] < gOflLumiCutoff &&
            this_channel_LB_current > gEpsilon) {
          pedestal_values.push_back(this_channel_LB_current);
        }
      }
      // Arbitrary restriction on how many good pedestal readings there must be
      const auto min_pedestal_values = 3;
      if (pedestal_values.size() < min_pedestal_values) {
        auto err_msg = "insufficient currents (< " +
                       std::to_string(min_pedestal_values) +
                       ") in channel " + channel_name +
                       " to form the avg pedestal calculation";
        return make_unexpected(make_shared<Error::Runtime>(err_msg, this_func_name));
      }

      auto pedestal_sum = std::accumulate(pedestal_values.begin(),
                                          pedestal_values.end(),
                                          0.0);
      auto pedestal_avg = pedestal_sum / pedestal_values.size();
      pedestals_.insert({{channel_name, pedestal_avg}});

      ++iChannel;
      pedestal_values.clear();
    }
  }

  if (analysis_->retrieve_currents()) {
    assert(currents_.size() == 0);

    auto iChannel = 0;
    auto this_channel_currents = vector<Float_t>();
    this_channel_currents.reserve(upper_LB_bound - lower_LB_bound + 1);
    for (const auto& channel: pedestals_) {
      auto this_channel_name = channel.first;
      auto this_channel_pedestal = channel.second;
      if (analysis_->use_baseline_subtraction_from_fit()) {
        this_channel_pedestal -= channel_calibrations_.at(this_channel_name).intercept/channel_calibrations_.at(this_channel_name).slope;
      }

      for (int i_LB = lower_LB_bound; i_LB <= upper_LB_bound; ++i_LB) {
        auto this_current = currents_temp[iChannel][i_LB] -
                            this_channel_pedestal;
        this_channel_currents.push_back(this_current);
      }

      currents_.insert({{std::move(this_channel_name),
                         std::move(this_channel_currents)}});
      ++iChannel;
      this_channel_currents.clear();
    }
  }

  if (analysis_->retrieve_lumi_ofl()) {
    assert(lumi_ofl_.size() == 0);

    // Mu-dependence correction evaluated by Vincent Hedberg and shown in this presentation by
    // Benedetto Giacobbe:
    // https://indico.cern.ch/event/435624/contributions/2209977/attachments/1293235/1927107/LTF_16062016.pdf
    auto apply_LUCID_mu_corr = analysis_->params().get<bool>("apply_LUCID_mu_corr");
    if (apply_LUCID_mu_corr) {
      auto conversion_factor = analysis_->x_sec() / (n_bunches_ * analysis_->f_rev());
      for (int i_LB = lower_LB_bound; i_LB <= upper_LB_bound; ++i_LB) {
        auto lumi = lumi_ofl_temp[i_LB];
        // have to convert to mu, apply the correction, then convert back to lumi
        auto mu = lumi*conversion_factor;
        auto mu_corr = -0.002*mu*mu + 1.008*mu; // this is the correction
        auto lumi_corr = mu_corr/conversion_factor;
        lumi_ofl_.push_back( lumi_corr );
      }
    }
    else {
      for (int i_LB = lower_LB_bound; i_LB <= upper_LB_bound; ++i_LB) {
        lumi_ofl_.push_back( lumi_ofl_temp[i_LB] );
      }
    }
  }

  return Void();
}

Expected<std::array<Int_t,2>> Run::GetLBBounds() const
{
  std::array<Int_t, 2> LB_bounds;
  const auto& custom_LB_bounds = analysis_->custom_LB_bounds();
  const auto target_LB_bounds = custom_LB_bounds.find(run_name_);
  if (target_LB_bounds != custom_LB_bounds.end()) {
    // LB number is 1-indexed, but we use 0-indexing to store it in memory
    // Negative value set by user is a placeholder for no desired bound to set
    // i.e. [-1,567] means use 567 as the upper bound and whatever the first RFP
    // LB is for the lower bound
    auto lower_bound =  target_LB_bounds->second.at(0) - 1;
    auto upper_bound =  target_LB_bounds->second.at(1) - 1;
    // i.e. unless we've specified a custom upper AND lower bound in config
    if (lower_bound < 0 || upper_bound < lower_bound) {
      auto RFP_LB_bounds = FindReadyForPhysicsLBBounds(RFP_flag_);
      RETURN_IF_ERR( RFP_LB_bounds )
      if (lower_bound < 0) {
        lower_bound = RFP_LB_bounds->at(0);
      }
      if (upper_bound < lower_bound) {
        upper_bound = RFP_LB_bounds->at(1);
      }
    }
    LB_bounds[0] = lower_bound;
    LB_bounds[1] = upper_bound;
  }
  else {
    auto RFP_LB_bounds = FindReadyForPhysicsLBBounds(RFP_flag_);
    RETURN_IF_ERR( RFP_LB_bounds )
    LB_bounds = *RFP_LB_bounds;
  }
  return LB_bounds;
}

Expected<Void> Run::CreateBenedettoOutput() const
{
  auto this_func_name = "Run::CreateBenedettoOutput";

  auto output_dir = analysis_->benedetto_output_dir();
  TRY( Util::mkdir(output_dir) )
  auto output_filepath = output_dir+run_name_+".dat";
  std::ofstream output_file(output_filepath);
  if (!output_file.is_open()) {
    return make_unexpected(make_shared<Error::File>(output_filepath, this_func_name));
  }
  cout << "        Writing LAr mu data to " << output_filepath << endl;

  auto num_points = mu_LAr_A_.size();
  output_file.precision(4);
  for (unsigned i_LB = 0; i_LB < num_points; ++i_LB) {
    output_file << (i_LB + LB_stability_offset_) << ' ' << mu_LAr_A_.at(i_LB) << ' '
                << mu_LAr_C_.at(i_LB) << std::endl;

    //auto mu_A = mu_LAr_A_.at(i_LB);
    //auto mu_C = mu_LAr_C_.at(i_LB);
    //Float_t conversion_factor = analysis_->x_sec() / (n_bunches_ * analysis_->f_rev());
    //auto mu_ofl = lumi_ofl_.at(i_LB) * conversion_factor;
    //output_file << (i_LB + LB_stability_offset_) << ' '
    //            << std::setw(5) << std::setfill('0') << std::left << std::noshowpos << mu_ofl << ' '
    //            << std::setw(5) << std::setfill('0') << std::left << mu_A << ' '
    //            << std::setw(5) << std::setfill('0') << std::left << mu_C << ' '
    //            << std::setw(5) << std::setfill('0') << std::left << std::showpos << 100*((mu_A + mu_C)/2 - mu_ofl) / mu_ofl << ' '
    //            << std::setw(5) << std::setfill('0') << std::left << 100*(mu_A - mu_C)/mu_A << std::endl;

    //output_file << (i_LB + LB_stability_offset_) << ' ' << lumi_ofl_.at(i_LB)/n_bunches_ << std::endl;
    //auto mu_A = mu_LAr_A_.at(i_LB);
    //auto mu_C = mu_LAr_C_.at(i_LB);
    //output_file << std::showpos << std::setw(3) << (i_LB + LB_stability_offset_) << "   " << 100*(mu_A - mu_C)/(mu_A + mu_C) << std::endl;
  }
  return Void();
}

void Run::GetExternalNBunches()
{
  auto& n_bunches = analysis_->n_bunches_from_file();
  auto res = n_bunches.find(run_name_);
  if (res != n_bunches.end()) {
    if (analysis_->verbose()) {
      if (n_bunches_) cout << "Warning: overriding non-negative n_bunches" << endl;
      cout << "Manually setting n_bunches for run " << run_name_ << endl;
    }
    n_bunches_ = res->second;
  }
  else {
    if (!n_bunches_) {
      if (analysis_->verbose()) {
        cout << "n_bunches = 0 for run " << run_name_ << 
          " and no value present in external file" << endl;
      }
    }
  }
}

// Calculates LAr luminosity from currents data for a run
Expected<Void> Run::CalcLArLumi()
{
  auto this_func_name = "Run::CalcLArLumi";
  //if (analysis_->verbose()) cout << "Calculating LAr luminosity" << endl;

  // Checks that these values were not previously calculated
  if (lumi_LAr_A_.size() > 0 || lumi_LAr_C_.size() > 0) {
    auto err_msg = "LAr lumi vector(s) already filled";
    return make_unexpected(make_shared<Error::Runtime>(err_msg, this_func_name));
  }

  // Checks that currents data exists
  if (currents_.size() == 0) {
    auto err_msg = "currents map empty";
    return make_unexpected(make_shared<Error::Runtime>(err_msg, this_func_name));
  }

  // All currents_ vectors should be the same size
  unsigned n_LB_stable = currents_.begin()->second.size();
  if (n_LB_stable == 0) {
    auto err_msg = "no stable lumi blocks";
    return make_unexpected(make_shared<Error::Runtime>(err_msg, this_func_name));
  }

  // Computes the beamspot correction based on given fit parameters and average beamspot z-position
  auto beamspot_corr_A = 1.0;
  auto beamspot_corr_C = 1.0;
  if (analysis_->use_beamspot_corr()) {
    auto intercept = analysis_->beamspot_corr_params().at(0);
    auto slope = analysis_->beamspot_corr_params().at(1);
    // beamspot correction parameters are in terms of % difference between A and C
    auto AC_percent_diff = intercept + slope*avg_beamspot_z_;
    beamspot_corr_C = AC_percent_diff/2/100 + 1;
    beamspot_corr_A = beamspot_corr_C/(AC_percent_diff/100 + 1);
    cout << beamspot_corr_C << '\t' << beamspot_corr_A << endl;
  }

  // Correction for luminosity-dependent effect seen in FCal (banana plot)
  auto apply_FCal_lumi_dep_corr = analysis_->params().get<bool>("apply_FCal_lumi_dep_corr");
  auto FCal_lumi_dep_corr = 1.0;

  // Averages lumi measurement from all channels for each side and for each
  //   lumi block
  for (unsigned i_LB = 0; i_LB < n_LB_stable; ++i_LB) {
    Float_t lumi_A_temp = 0.0;
    Float_t lumi_C_temp = 0.0;
    unsigned num_channels_A = 0;
    unsigned num_channels_C = 0;

    // Calculates lumi for a channel using its current and calibration and
    //   adds to running total
    for (const auto& this_channel_currents: currents_) {
      const auto& channel_name = this_channel_currents.first;
      auto current = this_channel_currents.second.at(i_LB);

      // Skip channel if current is ~0
      if (current < gLArCurrentCutoff) continue;

      auto intercept = analysis_->channel_calibrations().at(channel_name).intercept;
      auto slope = analysis_->channel_calibrations().at(channel_name).slope;

      if (FCalRegion::ToZSide(channel_name) == FCalRegion::ZSide::C) {
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
      lumi_LAr_A_.push_back(0.0);
    }
    else {
      auto lumi_LAr_A_this_LB = analysis_->anchoring_factor_A()*beamspot_corr_A*lumi_A_temp / num_channels_A;
      if (apply_FCal_lumi_dep_corr) {
        FCal_lumi_dep_corr = 0.99552064 + 0.12015841*std::exp(-1.04597906*lumi_LAr_A_this_LB/1000);
      }
      lumi_LAr_A_.push_back(lumi_LAr_A_this_LB/FCal_lumi_dep_corr);
    }
    if (num_channels_C == 0) {
      lumi_LAr_C_.push_back(0.0);
    }
    else {
      auto lumi_LAr_C_this_LB = analysis_->anchoring_factor_C()*beamspot_corr_C*lumi_C_temp / num_channels_C;
      if (apply_FCal_lumi_dep_corr) {
        FCal_lumi_dep_corr = 0.99552064 + 0.12015841*std::exp(-1.04597906*lumi_LAr_C_this_LB/1000);
      }
      lumi_LAr_C_.push_back(lumi_LAr_C_this_LB/FCal_lumi_dep_corr);
    }
  }

  return Void();
}

// Calculates <mu> from LAr luminosity for a run
Expected<Void> Run::CalcLArMu()
{
  auto this_func_name = "Run::CalcLArMu";

  //if (analysis_->verbose()) cout << "Calculating LAr <mu>" << endl;

  // Checks that LAr lumi data exists
  if (lumi_LAr_A_.size() == 0 && lumi_LAr_C_.size() == 0) {
    auto err_msg = "LAr lumi vector(s) have not been filled";
    return make_unexpected(make_shared<Error::Runtime>(err_msg, this_func_name));
  }
  // Checks that mu data has not already been calculated
  if (mu_LAr_A_.size() > 0 || mu_LAr_C_.size() > 0) {
    auto err_msg = "LAr mu vector(s) already filled";
    return make_unexpected(make_shared<Error::Runtime>(err_msg, this_func_name));
  }
  // Checks that value for number of collisions is nonzero
  if (n_bunches_ == 0) {
    auto err_msg = "num_collisions == 0";
    return make_unexpected(make_shared<Error::Runtime>(err_msg, this_func_name));
  }

  Float_t conversion_factor = analysis_->x_sec() / (n_bunches_ * analysis_->f_rev());

  // Calculates <mu> for each lumi block
  auto n_LB_stable = lumi_ofl_.size();
  mu_LAr_A_.reserve(n_LB_stable);
  for (auto lumi: lumi_LAr_A_) {
    mu_LAr_A_.push_back(lumi*conversion_factor);
  }
  mu_LAr_C_.reserve(n_LB_stable);
  for (auto lumi: lumi_LAr_C_) {
    mu_LAr_C_.push_back(lumi*conversion_factor);
  }
  mu_ofl_.reserve(n_LB_stable);
  for (auto lumi: lumi_ofl_) {
    mu_ofl_.push_back(lumi*conversion_factor);
  }
  return Void();
}

Expected<Void> Run::CreateLumiCurrentPlots() const
{
  auto this_func_name = "Run::CreateLumiCurrentPlots";
  if (analysis_->verbose()) cout << "    " << "Making lumi vs. current plots" << endl;

  LumiCurrentPlotOptions plot_options(analysis_->params(), "mode_options.lumi_current.");
  plot_options.run_name(string(run_name_));
  map<string, FitResults> fit_results;

  if (plot_options.do_individual()) {
    for (const auto &channel: currents_) {
      auto this_channel_name = channel.first;
      plot_options.channel_name(string(this_channel_name));

      auto channel_name = plot_options.channel_name();
      plot_options.title("Run "+plot_options.run_name()+", Channel "+channel_name);

      auto lumi_current_points = PointVectorFromVectors(
          lumi_ofl_,
          channel.second);

      CONTINUE_IF_ERR(lumi_current_points);

      auto points_filtered = *lumi_current_points;
      for (auto& point: points_filtered) {
        point[0] *= plot_options.x().scale();
        point[1] *= plot_options.y().scale();
      }

      points_filtered.erase(
          std::remove_if(
              points_filtered.begin(),
              points_filtered.end(),
              [] (auto point) {
                return point[0] < gOflLumiCutoff || point[1] < gLArCurrentCutoff;
              }
              ),
          points_filtered.end());

      auto this_channel_fit_results =
        Plotter::PlotLumiCurrent(points_filtered,
                                 plot_options);

      CONTINUE_IF_ERR(this_channel_fit_results)

      if (pedestals_.at(this_channel_name) > 20.0) {
        this_channel_fit_results->is_short = true;
      }
      else {
        this_channel_fit_results->is_short = false;
      }

      if (plot_options.do_fit()) {
        fit_results.insert({{this_channel_name, *this_channel_fit_results}});
      }
    } //channels loop
  }

  //if (plot_options.do_sum()) {
  //  vector<Float_t> channel_currents_sum_A;
  //  vector<Float_t> channel_currents_sum_C;
  //  for (const auto &channel: analysis_->channel_calibrations()) {
  //    auto this_channel_name = channel.first;
  //    auto this_channel_current = currents_.at(this_channel_name);

  //    vector<Float_t> *this_side_channel_currents_sum;
  //    if (this_channel_name.at(1) == '1') {
  //      this_side_channel_currents_sum = &channel_currents_sum_A;
  //    }
  //    else if (this_channel_name.at(1) == '8') {
  //      this_side_channel_currents_sum = &channel_currents_sum_C;
  //    }
  //    else {
  //      return make_unexpected(make_shared<Error::Logic>("Invalid channel name", this_func_name));
  //    }

  //    if (this_side_channel_currents_sum->size() == 0) {
  //      *this_side_channel_currents_sum = this_channel_current;
  //    }
  //    else {
  //      std::transform(this_side_channel_currents_sum->begin(),
  //                     this_side_channel_currents_sum->end(),
  //                     this_channel_current.begin(),
  //                     this_side_channel_currents_sum->begin(),
  //                     std::plus<Float_t>());
  //    }
  //  }

  //  // TODO: figure out how to deal with sums
  //  auto points_sum_A = PointVectorFromVectors(lumi_ofl_,
  //                                             channel_currents_sum_A);
  //  plot_options.channel_name("Sum_A");
  //  if (points_sum_A.valid()) {
  //    TRY( Plotter::PlotLumiCurrent(*points_sum_A, plot_options) )
  //  }

  //  auto points_sum_C = PointVectorFromVectors(lumi_ofl_,
  //                                             channel_currents_sum_C);
  //  plot_options.channel_name("Sum_C");
  //  if (points_sum_C.valid()) {
  //    TRY( Plotter::PlotLumiCurrent(*points_sum_C, plot_options) )
  //  }
  //}

  if (plot_options.do_fit()) {
    LOG_IF_ERR( Plotter::WriteFitResultsToTree(fit_results, plot_options) )
    LOG_IF_ERR( Plotter::WriteCalibrationToText(fit_results, plot_options) )
    //LOG_IF_ERR( Plotter::GeometricAnalysisOfFitResults(fit_results, plot_options) )
  }

  return Void();
}

// Creates those plots which use data from only a single sample
Expected<Void> Run::CreateSingleRunPlots() const
{
  for (const auto &mode: analysis_->modes()) {
    if (mode == "lumi_current") {
      if (analysis_->verbose()) {
        cout << "Creating plots for sample " << run_name_ << endl;
      }
      TRY( CreateLumiCurrentPlots() )
    }
  }
  // Remove the TRY above if other single run plots are added
  return Void();
}
