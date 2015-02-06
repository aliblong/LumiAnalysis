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
  TRY_THROW( PrepareAnalysis() )
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

Expected<Void> Analysis::PrepareAnalysis() {
// Reads in the analysis parameters (json) and channel information (text). 
// Sets flags for which data to retrieve based on which plots to produce.

  ReadParams();
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
    TRY( ReadCalibrations() )
  }

  return Void();
}

// Reads in calibration values for each of the channels being used (those
//   read in with ReadChannels)
Expected<Void> Analysis::ReadCalibrations() {
  auto this_func_name = "Analysis::ReadCalibrations";
  for (auto &channel: channel_calibrations_) {
    std::ifstream calibrations_file(calibrations_filepath_);
    if (!calibrations_file) {
      return make_unexpected(make_shared<Error::File>(calibrations_filepath_, this_func_name));
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

// Reads in parameters from json file and assigns their values to member
//   variables
void Analysis::ReadParams()
{
  JSONReader parameter_file(params_filepath_);

  verbose_ = parameter_file.get<bool>("verbose");

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

  calibrations_filepath_ = parameter_file.get<string>("input_filepaths.calibrations");
  channels_list_filepath_ =
    parameter_file.get<string>("input_filepaths.channels_list");
  pedestals_dir_ = parameter_file.get<string>("input_filepaths.pedestals");
  trees_dir_ = parameter_file.get<string>("input_filepaths.trees");
  run_list_dir_ = parameter_file.get<string>("input_filepaths.run_list");

  auto plot_types_map = parameter_file.get_map<bool>("plot_types");
  for (const auto &plot_type: plot_types_map) {
    if (plot_type.second == true) plot_types_.push_back(plot_type.first);
  }

  use_start_of_fill_pedestals_ = parameter_file.get<bool>("use_start_of_fill_pedestals");

  // Text output file for Benedetto
  do_benedetto_ = parameter_file.get<bool>("do_benedetto");
  benedetto_output_dir_ = parameter_file.get<string>("output_dirs.benedetto");

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

    TRY_CONTINUE( this_run.CreateSingleRunPlots() )

    if (retrieve_lumi_FCal_) {
      TRY_CONTINUE( this_run.CalcFCalLumi() )
      TRY_CONTINUE( this_run.CalcFCalMu() )
    }
    if (do_benedetto_) {
      TRY_CONTINUE( this_run.CreateBenedettoOutput() )
    }

    runs_data.insert({run_name, std::move(this_run)});
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
