#include "analysis.h"

#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include <Rtypes.h>

#include <boost/expected/expected.hpp>
#include <boost/container/flat_map.hpp>

#include "cutoffs.h"
#include "detector.h"
#include "json_reader.h"
#include "plotter.h"
#include "point.h"
#include "lumi_current_plot_options.h"
#include "mu_stab_plot_options.h"
#include "single_run_data.h"

using std::string;
using std::vector;
template<typename K, typename V>
using map = boost::container::flat_map<K, V>;

using std::cout;
using std::cerr;
using std::endl;

using std::make_shared;

using boost::expected;
using boost::make_unexpected;

using Error::Expected;

namespace {

Expected<Void> VerifyLBBounds(map<string, vector<int>> bounds)
{
  auto this_func_name = "VerifyLBBounds";
  for (const auto& item: bounds) {
    if (item.second.size() != 2) {
      auto err_msg = "Illegal number (!=2) of LB bounds for run " + item.first;
      return make_unexpected(make_shared<Error::Runtime>(err_msg, this_func_name));
    }
  }
  return Void();
}

Expected<vector<string>> VectorFromFile(const string& filepath)
{
  auto this_func_name = "VectorFromFile";

  std::ifstream file(filepath);
  if (!file) {
    return make_unexpected(make_shared<Error::File>(filepath, this_func_name));
  }

  string item;
  vector<string> items_vec;
  while ( getline(file, item) ) {
    // Skip samples commented out with #
    if (item[0] == '#') {
      cout << "Skipping " << item << endl;
      continue;
    }
    items_vec.push_back(item);
  }
  return items_vec;
}

  /*
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
*/

VectorP<Float_t> GenerateMuRatioVsMuPoints(const map<string, SingleRunData> &runs_data)
{
  VectorP<Float_t> points;
  // Roughly 500 points per run
  auto num_points_to_reserve = 500*2*runs_data.size();
  points.reserve(num_points_to_reserve);
  for (const auto& run: runs_data) {
    const auto& run_data = run.second;
    const auto& mu_A = run_data.mu_FCal_A();
    const auto& mu_C = run_data.mu_FCal_C();
    const auto& mu_ofl = run_data.mu_ofl();
    auto num_points = mu_ofl.size();
    for (auto i = 0; i < num_points; ++i) {
      auto mu_A_this_LB = mu_A[i];
      auto mu_C_this_LB = mu_C[i];
      auto mu_avg_this_LB = 0.0;
      if (mu_A_this_LB < gEpsilon) {
        mu_avg_this_LB = mu_C_this_LB;
      }
      else if (mu_C_this_LB < gEpsilon) {
        mu_avg_this_LB = mu_A_this_LB;
      }
      else {
        mu_avg_this_LB = (mu_A_this_LB + mu_C_this_LB)/2;
      }
      auto mu_ofl_this_LB = mu_ofl[i];
      auto mu_ratio_this_LB = (mu_avg_this_LB/mu_ofl_this_LB - 1)*100;
      points.push_back(Point<Float_t>{mu_ofl_this_LB, mu_ratio_this_LB});
    }
  }
  return points;
}

VectorP<Float_t> GenerateMuRatioVsLumiPoints(const map<string, SingleRunData> &runs_data)
{
  VectorP<Float_t> points;
  // Roughly 500 points per run
  auto num_points_to_reserve = 500*2*runs_data.size();
  points.reserve(num_points_to_reserve);
  for (const auto& run: runs_data) {
    const auto& run_data = run.second;
    const auto& lumi_A = run_data.lumi_FCal_A();
    const auto& lumi_C = run_data.lumi_FCal_C();
    const auto& lumi_ofl = run_data.lumi_ofl();
    auto num_points = lumi_ofl.size();
    for (auto i = 0; i < num_points; ++i) {
      auto lumi_A_this_LB = lumi_A[i];
      auto lumi_C_this_LB = lumi_C[i];
      auto lumi_FCal_avg_this_LB = 0.0;
      if (lumi_A_this_LB < gEpsilon) {
        lumi_FCal_avg_this_LB = lumi_C_this_LB;
      }
      else if (lumi_C_this_LB < gEpsilon) {
        lumi_FCal_avg_this_LB = lumi_A_this_LB;
      }
      else {
        lumi_FCal_avg_this_LB = (lumi_A_this_LB + lumi_C_this_LB)/2;
      }
      auto lumi_ofl_this_LB = lumi_ofl[i];
      auto lumi_ratio_this_LB = (lumi_FCal_avg_this_LB/lumi_ofl_this_LB - 1)*100;
      points.push_back(Point<Float_t>{lumi_ofl_this_LB, lumi_ratio_this_LB});
    }
  }
  return points;
}

VectorP<Float_t> GenerateLumiVsCurrentPoints(const map<string, SingleRunData> &runs_data)
{
  VectorP<Float_t> points;
  // Roughly 500 points per run
  auto num_points_to_reserve = 500*2*runs_data.size();
  points.reserve(num_points_to_reserve);
  for (const auto& run: runs_data) {
    const auto& run_data = run.second;
    const auto& lumi_A = run_data.lumi_FCal_A();
    const auto& lumi_C = run_data.lumi_FCal_C();
    const auto& lumi_ofl = run_data.lumi_ofl();
    auto num_points = lumi_ofl.size();
    for (auto i = 0; i < num_points; ++i) {
      auto lumi_A_this_LB = lumi_A[i];
      auto lumi_C_this_LB = lumi_C[i];
      auto lumi_FCal_avg_this_LB = 0.0;
      if (lumi_A_this_LB < gEpsilon) {
        lumi_FCal_avg_this_LB = lumi_C_this_LB;
      }
      else if (lumi_C_this_LB < gEpsilon) {
        lumi_FCal_avg_this_LB = lumi_A_this_LB;
      }
      else {
        lumi_FCal_avg_this_LB = (lumi_A_this_LB + lumi_C_this_LB)/2;
      }
      auto lumi_ofl_this_LB = lumi_ofl[i];
      points.push_back(Point<Float_t>{lumi_ofl_this_LB, lumi_FCal_avg_this_LB});
    }
  }
  return points;
}

} // Anonymous namespace

// Initializes various quantities required to run the analysis
Analysis::Analysis(string&& params_filepath)
{
  params_filepath_ = params_filepath;
  TRY_THROW( PrepareAnalysis() )
}

void Analysis::CreateAllRunPlots(const map<string, SingleRunData> &runs_data)
{
  for (const auto &plot_type: plot_types_) {
    if (plot_type == "mu_stability") {
      if (verbose_) cout << "Making mu stability plot" << endl;
      MuStabPlotOptions plot_options(params_filepath_);
      LOG_IF_ERR( Plotter::PlotMuStability(runs_data, plot_options) );
    }
    else if (plot_type == "mu_lumi_dependence") {
      if (verbose_) cout << "Making mu and lumi dependence plots" << endl;
      auto mu_dep_node = "plot_options.mu_dependence.";
      MuLumiDepPlotOptions mu_dep_plot_options(params_filepath_, mu_dep_node);
      auto mu_dep_points = GenerateMuRatioVsMuPoints(runs_data);
      LOG_IF_ERR( Plotter::PlotMuLumiDependence(mu_dep_points, mu_dep_plot_options) );

      auto lumi_dep_node = "plot_options.lumi_dependence.";
      MuLumiDepPlotOptions lumi_dep_plot_options(params_filepath_, lumi_dep_node);
      auto lumi_dep_points = GenerateMuRatioVsLumiPoints(runs_data);
      LOG_IF_ERR( Plotter::PlotMuLumiDependence(lumi_dep_points, lumi_dep_plot_options) );
    }
  }
}

// Reads in the analysis parameters (json) and channel information (text).
// Sets flags for which data to retrieve based on which plots to produce.
Expected<Void> Analysis::PrepareAnalysis()
{
  TRY( ReadParams() )
  for (const auto &plot_type: plot_types_) {
    if (plot_type == "lumi_current" ) {
      retrieve_currents_ = true;
      retrieve_lumi_ofl_ = true;
    }
    else if (plot_type == "mu_stability") {
      retrieve_timestamps_ = true;
      retrieve_lumi_ofl_ = true;
      retrieve_lumi_FCal_ = true;
    }
    else if (plot_type == "mu_lumi_dependence") {
      retrieve_lumi_ofl_ = true;
      retrieve_lumi_FCal_ = true;
    }
  }

  if (do_benedetto_) {
    retrieve_lumi_FCal_ = true;
    //TODO: remove this
    retrieve_lumi_ofl_ = true;
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
Expected<Void> Analysis::ReadCalibrations()
{
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
Expected<Void> Analysis::ReadChannels()
{
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
    channel_calibrations_.emplace_hint(end(channel_calibrations_),
                          std::move(channel_name), ChannelCalibration{0.0, 0.0});
  }
  channels_list_file.close();

  return Void();
}

// Reads in parameters from json file and assigns their values to member
//   variables
Error::Expected<Void> Analysis::ReadParams()
{
  JSONReader parameter_file(params_filepath_);

  auto detector = Detector::FromString(parameter_file.get<string>("detector"));
  RETURN_IF_ERR( detector )
  detector_ = *detector;

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

  JSONReader LB_bounds_file(parameter_file.get<string>("input_filepaths.custom_LB_bounds"));
  auto custom_LB_bounds = LB_bounds_file.get_map_of_vectors<int>("");

  TRY( VerifyLBBounds(custom_LB_bounds) )
  custom_LB_bounds_ = custom_LB_bounds;

  // Text output file for Benedetto
  do_benedetto_ = parameter_file.get<bool>("do_benedetto");
  benedetto_output_dir_ = parameter_file.get<string>("output_dirs.base") +
                          parameter_file.get<string>("output_dirs.benedetto");

  return Void();
}

// Control flow for the analysis of samples and creation of plots
Expected<Void> Analysis::RunAnalysis()
{
  auto run_names_vec = VectorFromFile(run_list_dir_);
  RETURN_IF_ERR( run_names_vec )

  // Iterate over samples
  map<string, SingleRunData> runs_data;

  for (const auto &run_name: *run_names_vec) {
    auto this_run = SingleRunData(run_name, this);
    TRY_CONTINUE( this_run.Init() )

    TRY_CONTINUE( this_run.CreateSingleRunPlots() )
    if (retrieve_lumi_FCal_) {
      TRY_CONTINUE( this_run.CalcFCalLumi() )
      TRY_CONTINUE( this_run.CalcFCalMu() )
    }
    if (do_benedetto_) {
      TRY_CONTINUE( this_run.CreateBenedettoOutput() )
    }

    runs_data.emplace_hint(end(runs_data),  //change this to cend in C++14
                           string(run_name), std::move(this_run));
  }

  CreateAllRunPlots(runs_data);

  return Void();
}
