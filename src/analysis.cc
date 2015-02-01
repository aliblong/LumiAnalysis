#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "Rtypes.h"

#include "boost/expected/expected.hpp"

#include "analysis.h"
#include "cutoffs.h"
#include "json_reader.h"
#include "plotter.h"
#include "point.h"
#include "lumi_current_plot_options.h"
#include "mu_stab_plot_options.h"
#include "single_run_data.h"

using std::string;
using std::vector;
using std::map;

using std::cout;
using std::cerr;
using std::endl;

using std::make_shared;

using boost::expected;
using boost::make_unexpected;

using Error::Expected;

namespace {

void DumpVector(const vector<Float_t> &vec) {
  // Prints newline-separated vector elements to cout.
  for (auto element: vec) cout << element << endl;
}

// Checks whether or not the vector has all elements less than epsilon
bool IsZeroes(const vector<Float_t> &vec) {
  for (auto element: vec) {
    if (element > gEpsilon) return false;
  }
  return true;
}

} // Anonymous namespace

// Initializes various quantities required to run the analysis
Analysis::Analysis(string params_filepath)
  : verbose_(false),
    do_benedetto_(false),

    f_rev_(0.0),
    x_sec_(0.0),
    ref_run_number_(0),
    corr_A_(0.0),
    corr_C_(0.0),
    corr_Avg_(0.0),

    params_filepath_(""),
    calibrations_filepath_(""),
    channels_list_filepath_(""),
    pedestals_dir_(""),
    trees_dir_(""),
    run_list_dir_(""),
    retrieve_timestamps_(false),
    retrieve_currents_(false),
    retrieve_lumi_BCM_(false),
    retrieve_lumi_FCal_(false) {

  params_filepath_ = params_filepath;
  TRY_THROW( PrepareAnalysis(params_filepath) )
}

// Calculates FCal luminosity from currents data for a run
Expected<Void> Analysis::CalcFCalLumi(SingleRunData &this_run) {
  auto this_func_name = "Analysis::CalcFCalLumi";
  if (verbose_) cout << "Calculating FCal luminosity" << endl;

  // Checks that these values were not previously calculated
  if (this_run.lumi_FCal_A_.size() > 0 || this_run.lumi_FCal_C_.size() > 0) {
    auto err_msg = "FCal lumi vector(s) already filled";
    return make_unexpected(make_shared<Error::Runtime>(err_msg, this_func_name));
  }

  // Checks that currents data exists
  if (this_run.currents_.size() == 0) {
    auto err_msg = "currents map empty";
    return make_unexpected(make_shared<Error::Runtime>(err_msg, this_func_name));
  }

  unsigned nLB_stable = this_run.currents_.begin()->second.size();
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
    for (const auto &this_module_currents: this_run.currents_) {
      auto current = this_module_currents.second.at(iLB);

      // Skip channel if current is too low
      if (current < gFCalCurrentCutoff) continue;

      string module_name = this_module_currents.first;
      auto intercept = channel_calibrations_.at(module_name).intercept;
      auto slope = channel_calibrations_.at(module_name).slope;

      if (module_name.at(1) == '1') { // module from C-side
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
      this_run.lumi_FCal_A_.push_back(0.0);
    }
    else {
      this_run.lumi_FCal_A_.push_back(corr_A_*lumi_A_temp / num_channels_A);
    }
    if (num_channels_C == 0) {
      this_run.lumi_FCal_C_.push_back(0.0);
    }
    else {
      this_run.lumi_FCal_C_.push_back(corr_C_*lumi_C_temp / num_channels_C);
    }
  }

  return Void();
}

// Calculates <mu> from FCal luminosity for a run
Expected<Void> Analysis::CalcFCalMu(SingleRunData &this_run) {
  auto this_func_name = "Analysis::CalcFCalMu";

  if (verbose_) cout << "Calculating FCal <mu>" << endl;

  // Checks that FCal lumi data exists
  if (this_run.lumi_FCal_A_.size() == 0 && this_run.lumi_FCal_C_.size() == 0) {
    auto err_msg = "FCal lumi vector(s) have not been filled";
    return make_unexpected(make_shared<Error::Runtime>(err_msg, this_func_name));
  }
  // Checks that mu data has not already been calculated
  if (this_run.mu_FCal_A_.size() > 0 || this_run.mu_FCal_C_.size() > 0) {
    auto err_msg = "FCal mu vector(s) already filled";
    return make_unexpected(make_shared<Error::Runtime>(err_msg, this_func_name));
  }
  // Checks that value for number of collisions is nonzero
  if (this_run.nCollisions_ == 0) {
    auto err_msg = "num_collisions == 0";
    return make_unexpected(make_shared<Error::Runtime>(err_msg, this_func_name));
  }

  Float_t conversion_factor = x_sec_ / (this_run.nCollisions_ * f_rev_);
  // Cross-section conversion factor from 7 -> 8 TeV
  if (this_run.run_name_.at(0) == '2') conversion_factor /= 1.05;

  // Calculates <mu> for each lumi block
  for (const auto &lumi: this_run.lumi_FCal_A_) {
    this_run.mu_FCal_A_.push_back(lumi*conversion_factor);
  }
  for (const auto &lumi: this_run.lumi_FCal_C_) {
    this_run.mu_FCal_C_.push_back(lumi*conversion_factor);
  }
  return Void();
}

void Analysis::CreateAllRunPlots(const map<string, SingleRunData> &runs_data) {
  for (const auto &plot_type: plot_types_) {
    if (plot_type == "mu_stability") {
      if (verbose_) cout << "Making mu stability plot" << endl;
      MuStabPlotOptions plot_options(params_filepath_);
      LOG_IF_ERR( Plotter::PlotMuStability(runs_data, plot_options) );
    }
  }
}

Expected<Void> Analysis::CreateLumiCurrentPlots(const SingleRunData &this_run) {
  const char* this_func_name = "Analysis::CreateLumiCurrentPlots";
  if (verbose_) cout << "    " << "Making lumi vs. current plots" << endl;

  LumiCurrentPlotOptions plot_options(params_filepath_);
  plot_options.set_run_name(this_run.run_name_);
  map<string, FitResults> fit_results;

  if (plot_options.do_individual()) {
    for (const auto &channel: channel_calibrations_) {
      auto this_channel_name = channel.first;
      plot_options.set_channel_name(this_channel_name);

      auto lumi_current_points = PointVectorFromVectors(
          this_run.lumi_BCM(),
          this_run.currents().at(this_channel_name));

      CONTINUE_IF_ERR(lumi_current_points);

      auto this_channel_fit_results =
        Plotter::PlotLumiCurrent(*lumi_current_points,
                                 plot_options);

      CONTINUE_IF_ERR(this_channel_fit_results)

      if (this_run.pedestals_.at(this_channel_name) > 20) {
        this_channel_fit_results->is_short = true;
      }
      else {
        this_channel_fit_results->is_short = false;
      }

      if (plot_options.do_fit()) {
        fit_results.insert(std::make_pair(this_channel_name,
                                          *this_channel_fit_results));
      }
    } //channels loop
  }

  if (plot_options.do_sum()) {
    vector<Float_t> channel_currents_sum_A;
    vector<Float_t> channel_currents_sum_C;
    for (const auto &channel: channel_calibrations_) {
      auto this_channel_name = channel.first;
      auto this_channel_current = this_run.currents_.at(this_channel_name);

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
    auto points_sum_A = PointVectorFromVectors(this_run.lumi_BCM(),
                                               channel_currents_sum_A);
    plot_options.set_channel_name("Sum_A");
    if (points_sum_A.valid()) {
      TRY( Plotter::PlotLumiCurrent(*points_sum_A, plot_options) )
    }

    auto points_sum_C = PointVectorFromVectors(this_run.lumi_BCM(),
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
Expected<Void> Analysis::CreateSingleRunPlots(const SingleRunData &this_run) {
  string run_name = this_run.run_name_;

  for (const auto &plot_type: plot_types_) {
    if (plot_type == "lumi_current") {
      if (verbose_) cout << "Creating plots for sample " << run_name << endl;
      map<string, FitResults> fit_results;
      TRY( CreateLumiCurrentPlots(this_run) )
    }
  }
  // Remove the TRY above if other single run plots are added
  return Void();
}

Expected<Void> Analysis::PrepareAnalysis(string params_filepath) {
// Reads in the analysis parameters (json) and channel information (text). 
// Sets flags for which data to retrieve based on which plots to produce.

  ReadParams(params_filepath);
  for (const auto &plot_type: plot_types_) {
    if (plot_type == "lumi_current" ) {
      retrieve_currents_ = true;
      retrieve_lumi_BCM_ = true;
    }
    else if (plot_type == "mu_stability") {
      retrieve_timestamps_ = true;
      retrieve_lumi_BCM_ = true;
      retrieve_lumi_FCal_ = true;
    }
  }

  if (do_benedetto_) {
    retrieve_lumi_FCal_ = true;
  }

  // FCal currents are required to calculate FCal lumi
  if (retrieve_lumi_FCal_) retrieve_currents_ = true;

  TRY( ReadChannels() )

  if (retrieve_lumi_FCal_) {
    TRY( ReadCalibrations(calibrations_filepath_) )
  }

  return Void();
}

// Reads in calibration values for each of the channels being used (those
//   read in with ReadChannels)
Expected<Void> Analysis::ReadCalibrations(string calibrations_filepath) {
  auto this_func_name = "Analysis::ReadCalibrations";
  for (auto &channel: channel_calibrations_) {
    std::ifstream calibrations_file(calibrations_filepath);
    if (!calibrations_file) {
      return make_unexpected(make_shared<Error::File>(calibrations_filepath, this_func_name));
    }

    string channel_name;
    Float_t slope;
    Float_t intercept;
    bool found_channel = false;
    while (calibrations_file >> channel_name >> slope >> intercept) {
      if (channel.first == channel_name) {
        channel.second.slope = slope;
        channel.second.intercept = intercept;
        found_channel = true;
        break;
      }
    }

    calibrations_file.close();

    if (!found_channel) {
      auto err_msg = "could not locate calibration for channel "+channel.first;
      return make_unexpected(make_shared<Error::Runtime>(err_msg, this_func_name));
    }
  }

  return Void();
}

// Reads the list of channels to be used in the analysis, and stores them in a
//   map
Expected<Void> Analysis::ReadChannels() {
  auto this_func_name = "Analysis::ReadChannels";

  std::ifstream channels_list_file(channels_list_filepath_);
  if (!channels_list_file) {
    return make_unexpected(make_shared<Error::File>(channels_list_filepath_, this_func_name));
  }

  string channel_name;
  while ( getline(channels_list_file, channel_name) ) {
    if (channel_name[0] == '#') {
      if (verbose_) cout << "Skipping channel " << channel_name << endl;
      continue;
    }
    channel_calibrations_.insert({channel_name, {0.0, 0.0}});
  }
  channels_list_file.close();

  return Void();
}

void Analysis::ReadParams(string params_filepath) {
// Reads in parameters from json file and assigns their values to member
//   variables

  JSONReader parameter_file(params_filepath);

  verbose_ = parameter_file.get<bool>("verbose");

  // Text output file for Benedetto
  do_benedetto_ = parameter_file.get<bool>("do_benedetto");
  benedetto_output_dir_ = parameter_file.get<string>("output_dirs.benedetto");

  // Used to calculate FCal luminosity
  //   Bunch crossing frequency
  f_rev_ = parameter_file.get<double>("lumi_calculation.f_rev");
  //   Cross-section for pp interaction @ 7 TeV
  x_sec_ = parameter_file.get<double>("lumi_calculation.x_sec");

  // Value of 0 for reference run means set no corrections
  ref_run_number_ = parameter_file.get<int>("ref_run_number");
  string ref_run_str = std::to_string(ref_run_number_);
  corr_A_ = parameter_file.get<double>("corrections."+ref_run_str+".A");
  corr_C_ = parameter_file.get<double>("corrections."+ref_run_str+".C");
  corr_Avg_ = parameter_file.get<double>("corrections."+ref_run_str+".Avg");

  auto plot_types_map = parameter_file.get_map<bool>("plot_types");
  for (const auto &plot_type: plot_types_map) {
    if (plot_type.second == true) plot_types_.push_back(plot_type.first);
  }

  calibrations_filepath_ = parameter_file.get<string>("input_filepaths.calibrations");
  channels_list_filepath_ =
    parameter_file.get<string>("input_filepaths.channels_list");
  pedestals_dir_ = parameter_file.get<string>("input_filepaths.pedestals");
  trees_dir_ = parameter_file.get<string>("input_filepaths.trees");
  run_list_dir_ = parameter_file.get<string>("input_filepaths.run_list");
}

// Control flow for the analysis of samples and creation of plots
Expected<Void> Analysis::RunAnalysis() {

  auto this_func_name = "Analysis::RunAnalysis()";

  std::ifstream run_names_file(run_list_dir_);
  if (!run_names_file) {
    return make_unexpected(make_shared<Error::File>(run_list_dir_, this_func_name));
  }

  string run_name;
  vector<string> run_names_vec;
  while ( getline(run_names_file, run_name) ) {
    // Skip samples commented out with #
    if (run_name[0] == '#') {
      if (verbose_) cout << "Skipping run " << run_name << endl;
      continue;
    }
    run_names_vec.push_back(run_name);
  }

  // Iterate over samples
  map<string, SingleRunData> runs_data;

  for (const auto &run_name: run_names_vec) {
    auto this_run = SingleRunData{run_name, *this};

    TRY_CONTINUE( CreateSingleRunPlots(this_run) )

    if (retrieve_lumi_FCal_) {
      TRY_CONTINUE( CalcFCalLumi(this_run) )
      TRY_CONTINUE( CalcFCalMu(this_run) )
    }
    if (do_benedetto_) {
      TRY_CONTINUE( this_run.CreateBenedettoOutput() )
    }

    runs_data.insert(std::make_pair(run_name, std::move(this_run)));
  }

  CreateAllRunPlots(runs_data);

  return Void();
}

/*
int WriteCurrentsToFile(vector<Float_t> currents, string run_name) {
  std::ofstream out_file(out_filepath);
  if (!out_file.is_open()) {
    cerr << "ERROR: Could not open file \'" << out_filepath << "\'" << endl;
    return 1;
  }
  for (auto iEvent: currents_) {
    for (auto iCurrent: iEvent) out_file << iCurrent << ' ';
    out_file << '\n';
  }
  out_file.close();
  return 0;
}

template <typename T>
int WriteVectorToFile(const vector<T> &vec, string filepath) {
  // Writes each element of a vector to a file at the specified path;
  //   elements are separated by newlines.

  std::ofstream out_file(filepath);
  if (!out_file.is_open()) {
    cerr << "ERROR: could not open file \'" << filepath << "\'" << endl;
    return 1;
  }

  for (auto iElement: vec) out_file << iElement << '\n';
  out_file.close();

  return 0;
}

*/
