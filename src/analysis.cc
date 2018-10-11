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
#include "run.h"
#include "util.h"

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

VectorP<Float_t> GenerateMuRatioVsMuPoints(
    const map<string, Run> &runs_data,
    double x_scale,
    double y_scale)
{
  VectorP<Float_t> points;
  // Roughly 500 points per run
  auto num_points_to_reserve = 500*2*runs_data.size();
  points.reserve(num_points_to_reserve);
  for (const auto& run: runs_data) {
    const auto& run_data = run.second;
    const auto& mu_A = run_data.mu_LAr_A();
    const auto& mu_C = run_data.mu_LAr_C();
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

VectorP<Float_t> GenerateMuRatioVsLumiPoints(
    const map<string, Run> &runs_data,
    double x_scale,
    double y_scale)
{
  VectorP<Float_t> points;
  // Roughly 500 points per run
  auto num_points_to_reserve = 500*2*runs_data.size();
  points.reserve(num_points_to_reserve);
  for (const auto& run: runs_data) {
    const auto& run_data = run.second;
    const auto& lumi_A = run_data.lumi_LAr_A();
    const auto& lumi_C = run_data.lumi_LAr_C();
    const auto& lumi_ofl = run_data.lumi_ofl();
    auto num_points = lumi_ofl.size();
    for (auto i = 0; i < num_points; ++i) {
      auto lumi_A_this_LB = lumi_A[i];
      auto lumi_C_this_LB = lumi_C[i];
      auto lumi_LAr_avg_this_LB = 0.0;
      if (lumi_A_this_LB < gEpsilon) {
        lumi_LAr_avg_this_LB = lumi_C_this_LB;
      }
      else if (lumi_C_this_LB < gEpsilon) {
        lumi_LAr_avg_this_LB = lumi_A_this_LB;
      }
      else {
        lumi_LAr_avg_this_LB = (lumi_A_this_LB + lumi_C_this_LB)/2;
      }
      auto lumi_ofl_this_LB = lumi_ofl[i];
      auto lumi_ratio_this_LB = (lumi_LAr_avg_this_LB/lumi_ofl_this_LB - 1)*100;
      points.push_back(Point<Float_t>{lumi_ofl_this_LB, lumi_ratio_this_LB});
    }
  }
  return points;
}

VectorP<Float_t> GenerateAvgMuRatioVsLumiPoints(
    const map<string, Run> &runs_data,
    double x_scale,
    double y_scale)
{
  VectorP<Float_t> points;
  // Roughly 500 points per run
  auto num_points_to_reserve = runs_data.size();
  points.reserve(num_points_to_reserve);
  for (const auto& run: runs_data) {
    const auto& run_data = run.second;
    const auto& lumi_A = run_data.lumi_LAr_A();
    const auto& lumi_C = run_data.lumi_LAr_C();
    const auto& lumi_ofl = run_data.lumi_ofl();
    auto num_points = lumi_ofl.size();
    auto lumi_ofl_sum_this_run = 0.0;
    auto lumi_ratio_sum_this_run = 0.0;
    auto n_LB_with_nonzero_lumi = 0;
    for (auto i = 0; i < num_points; ++i) {
      auto lumi_A_this_LB = lumi_A[i];
      auto lumi_C_this_LB = lumi_C[i];
      auto lumi_LAr_avg_this_LB = 0.0;
      auto lumi_ofl_this_LB = lumi_ofl[i];
      if (lumi_ofl_this_LB < gEpsilon) continue;
      if (lumi_A_this_LB < gEpsilon && lumi_C_this_LB < gEpsilon) {
        continue;
      }
      else if (lumi_A_this_LB < gEpsilon) {
        lumi_LAr_avg_this_LB = lumi_C_this_LB;
      }
      else if (lumi_C_this_LB < gEpsilon) {
        lumi_LAr_avg_this_LB = lumi_A_this_LB;
      }
      else {
        lumi_LAr_avg_this_LB = (lumi_A_this_LB + lumi_C_this_LB)/2;
      }
      auto lumi_ratio_this_LB = (lumi_LAr_avg_this_LB/lumi_ofl_this_LB - 1)*100;
      lumi_ofl_sum_this_run += lumi_ofl_this_LB;
      lumi_ratio_sum_this_run += lumi_ratio_this_LB;
      ++n_LB_with_nonzero_lumi;
    }
    if (n_LB_with_nonzero_lumi == 0) continue;
    auto avg_lumi_ofl_this_run = lumi_ofl_sum_this_run / n_LB_with_nonzero_lumi;
    auto avg_lumi_ratio_this_run = lumi_ratio_sum_this_run / n_LB_with_nonzero_lumi;
    points.push_back({avg_lumi_ofl_this_run*x_scale, avg_lumi_ratio_this_run*y_scale});
  }
  return points;
}

VectorP<Float_t> GenerateAvgMuRatioVsBeamspotZPoints(
    const map<string, Run> &runs_data,
    double x_scale,
    double y_scale)
{
  VectorP<Float_t> points;
  // Roughly 500 points per run
  auto num_points_to_reserve = runs_data.size();
  points.reserve(num_points_to_reserve);
  for (const auto& run: runs_data) {
    const auto& run_data = run.second;
    const auto& lumi_A = run_data.lumi_LAr_A();
    const auto& lumi_C = run_data.lumi_LAr_C();
    const auto& lumi_ofl = run_data.lumi_ofl();
    auto num_points = lumi_ofl.size();
    auto lumi_ofl_sum_this_run = 0.0;
    auto lumi_ratio_sum_this_run = 0.0;
    auto n_LB_with_nonzero_lumi = 0;
    for (auto i = 0; i < num_points; ++i) {
      auto lumi_A_this_LB = lumi_A[i];
      auto lumi_C_this_LB = lumi_C[i];
      auto lumi_LAr_avg_this_LB = 0.0;
      auto lumi_ofl_this_LB = lumi_ofl[i];
      if (lumi_ofl_this_LB < gEpsilon) continue;
      if (lumi_A_this_LB < gEpsilon && lumi_C_this_LB < gEpsilon) {
        continue;
      }
      else if (lumi_A_this_LB < gEpsilon) {
        lumi_LAr_avg_this_LB = lumi_C_this_LB;
      }
      else if (lumi_C_this_LB < gEpsilon) {
        lumi_LAr_avg_this_LB = lumi_A_this_LB;
      }
      else {
        lumi_LAr_avg_this_LB = (lumi_A_this_LB + lumi_C_this_LB)/2;
      }
      auto lumi_ratio_this_LB = (lumi_LAr_avg_this_LB/lumi_ofl_this_LB - 1)*100;
      lumi_ofl_sum_this_run += lumi_ofl_this_LB;
      lumi_ratio_sum_this_run += lumi_ratio_this_LB;
      ++n_LB_with_nonzero_lumi;
    }
    if (n_LB_with_nonzero_lumi == 0) continue;
    auto avg_lumi_ofl_this_run = lumi_ofl_sum_this_run / n_LB_with_nonzero_lumi;
    auto avg_lumi_ratio_this_run = lumi_ratio_sum_this_run / n_LB_with_nonzero_lumi;
    points.push_back({run_data.avg_beamspot_z()*x_scale, avg_lumi_ratio_this_run*y_scale});
  }
  return points;
}

VectorP<Float_t> GenerateACRatioVsBeamspotZPoints(
    const map<string, Run> &runs_data,
    double x_scale,
    double y_scale)
{
  VectorP<Float_t> points;
  auto num_points_to_reserve = runs_data.size();
  points.reserve(num_points_to_reserve);
  for (const auto& run: runs_data) {
    const auto& run_data = run.second;
    const auto& lumi_A = run_data.lumi_LAr_A();
    const auto& lumi_C = run_data.lumi_LAr_C();
    auto num_points = lumi_A.size();
    auto AC_ratio_sum = 0.0;
    auto n_LB_with_nonzero_A_and_C_lumi = 0;
    for (auto i = 0; i < num_points; ++i) {
      auto lumi_A_this_LB = lumi_A[i];
      auto lumi_C_this_LB = lumi_C[i];
      //cout << lumi_A_this_LB << '\t' << lumi_C_this_LB << endl;
      if (lumi_A_this_LB < gEpsilon || lumi_C_this_LB < gEpsilon) {
        continue;
      }
      AC_ratio_sum +=  lumi_A_this_LB / lumi_C_this_LB;
      ++n_LB_with_nonzero_A_and_C_lumi;
    }
    //cout << n_LB_with_nonzero_A_and_C_lumi << endl;
    if (n_LB_with_nonzero_A_and_C_lumi == 0) continue;
    auto avg_AC_ratio_percent_diff = (AC_ratio_sum / n_LB_with_nonzero_A_and_C_lumi - 1)*100;
    cout << run_data.avg_beamspot_z() << '\t' << avg_AC_ratio_percent_diff << endl;
    points.push_back({run_data.avg_beamspot_z()*x_scale, avg_AC_ratio_percent_diff*y_scale});
  }
  return points;
}

VectorP<Float_t> GenerateLumiVsCurrentPoints(
    const map<string, Run> &runs_data,
    const string& channel_name)
{
  VectorP<Float_t> points;
  // Roughly 500 points per run
  auto num_points_to_reserve = 500*2*runs_data.size();
  points.reserve(num_points_to_reserve);
  for (const auto& run: runs_data) {
    const auto& run_data = run.second;
    const auto& current = run_data.currents().at(channel_name);
    const auto& lumi_ofl = run_data.lumi_ofl();
    auto num_points = lumi_ofl.size();
    for (auto i = 0; i < num_points; ++i) {
      points.push_back(Point<Float_t>{lumi_ofl[i], current[i]});
    }
  }
  return points;
}

} // Anonymous namespace

// Initializes various quantities required to run the analysis
Analysis::Analysis(string&& params_filepath) :
  params_(params_filepath)
{
  params_.SetDefaults();
  TRY_THROW( PrepareAnalysis() )
}

void Analysis::CreateAllRunPlots(const map<string, Run> &runs_data)
{
  for (const auto &mode: modes()) {
    if (mode == "mu_stability") {
      VPRINT("Making mu stability plot")
      auto node = "mode_options.mu_stability.";
      MuStabPlotOptions plot_options(params_, node);
      LOG_IF_ERR( Plotter::PlotMuStability(runs_data, plot_options) );
    }
    else if (mode == "mu_lumi_dependence") {
      VPRINT("Making mu and lumi dependence plots")
      {
        auto node = "mode_options.mu_dependence.";
        MuLumiDepPlotOptions plot_options(params_, node);
        auto points = GenerateMuRatioVsMuPoints(
            runs_data,
            plot_options.x().scale(),
            plot_options.y().scale());
        LOG_IF_ERR( Plotter::PlotMuLumiDependence(points, plot_options) );
      }
      {
        auto node = "mode_options.lumi_dependence.";
        MuLumiDepPlotOptions plot_options(params_, node);
        auto points = GenerateMuRatioVsLumiPoints(
            runs_data,
            plot_options.x().scale(),
            plot_options.y().scale());
        LOG_IF_ERR( Plotter::PlotMuLumiDependence(points, plot_options) );
      }
      {
        auto node = "mode_options.avg_lumi_dependence.";
        MuLumiDepPlotOptions plot_options(params_, node);
        auto points = GenerateAvgMuRatioVsLumiPoints(
            runs_data,
            plot_options.x().scale(),
            plot_options.y().scale());
        LOG_IF_ERR( Plotter::PlotMuLumiDependence(points, plot_options) );
      }
    }
    else if (mode == "lumi_current_multirun") {
      VPRINT("Making current vs. luminosity plot for multiple runs")
      auto node = "mode_options.lumi_current.";
      LumiCurrentPlotOptions plot_options(params_, node);
      plot_options.run_name("multirun");
      map<string, FitResults> fit_results;
      for (const auto& channel: channel_calibrations()) {
        auto this_channel_name = channel.first;
        auto lumi_current_points = GenerateLumiVsCurrentPoints(runs_data, this_channel_name);
        auto is_verbose = verbose();
        lumi_current_points.erase(
            std::remove_if(
                lumi_current_points.begin(),
                lumi_current_points.end(),
                [is_verbose] (auto point) {
                  auto lumi = point[0];
                  auto current = point[1];
                  if (is_verbose && current < 0) {
                    // Placeholder current representing an absent value is -9999.0
                    if (fabs(current + 9999.0) < 1.0) {
                      cout << "No current for this LB" << endl;
                    }
                    else {
                      cout << "Negative current. Current == " << current << ", lumi == " << lumi << endl;
                    }
                  }
                  return lumi < gOflLumiCutoff || current < gLArCurrentCutoff;
                }
                ),
            lumi_current_points.end());
        plot_options.channel_name(std::move(this_channel_name));
        plot_options.title("Runs "+runs_data.begin()->first+" - "+runs_data.rbegin()->first+", Channel "+this_channel_name);
        auto this_channel_fit_results = Plotter::PlotLumiCurrent(lumi_current_points, plot_options);
        CONTINUE_IF_ERR(this_channel_fit_results)

        if (plot_options.do_fit()) {
          fit_results.insert({{this_channel_name, *this_channel_fit_results}});
        }
      }
      LOG_IF_ERR( Plotter::WriteFitResultsToTree(fit_results, plot_options) )
      LOG_IF_ERR( Plotter::WriteCalibrationToText(fit_results, plot_options) )
    }
    else if (mode == "beamspot_AC") {
      VPRINT("Making beamspot Z-position vs. A/C ratio plot")
      auto node = "mode_options.beamspot_AC.";
      MuLumiDepPlotOptions plot_options(params_, node);
      auto points = GenerateACRatioVsBeamspotZPoints(
          runs_data,
          plot_options.x().scale(),
          plot_options.y().scale());
      LOG_IF_ERR( Plotter::PlotMuLumiDependence(points, plot_options) );
    }
    else if (mode == "beamspot_LAr-ofl") {
      VPRINT("Making beamspot Z-position vs. LAr/ofl ratio plot")
      auto node = "mode_options.beamspot_LAr-ofl.";
      MuLumiDepPlotOptions plot_options(params_, node);
      auto points = GenerateAvgMuRatioVsBeamspotZPoints(
          runs_data,
          plot_options.x().scale(),
          plot_options.y().scale());
      LOG_IF_ERR( Plotter::PlotMuLumiDependence(points, plot_options) );
    }
  }
}

// Reads in the analysis parameters (json) and channel information (text).
// Sets flags for which data to retrieve based on which plots to produce.
Expected<Void> Analysis::PrepareAnalysis()
{
  for (const auto &mode: modes()) {
    if (mode == "lumi_current" || mode == "lumi_current_multirun") {
      retrieve_currents_ = true;
      retrieve_lumi_ofl_ = true;
    }
    else if (mode == "mu_stability") {
      retrieve_timestamps_ = true;
      retrieve_lumi_ofl_ = true;
      retrieve_lumi_LAr_ = true;
    }
    else if (mode == "mu_lumi_dependence") {
      retrieve_lumi_ofl_ = true;
      retrieve_lumi_LAr_ = true;
      retrieve_mu_LAr_ = true;
    }
    else if (mode == "beamspot_AC") {
      retrieve_lumi_LAr_ = true;
      retrieve_beamspot_ = true;
    }
    else if (mode == "beamspot_LAr-ofl") {
      retrieve_lumi_ofl_ = true;
      retrieve_lumi_LAr_ = true;
      retrieve_beamspot_ = true;
    }
    else if (mode == "benedetto") {
      retrieve_lumi_LAr_ = true;
      retrieve_mu_LAr_ = true;
      //TODO: remove this
      retrieve_lumi_ofl_ = true;
    }
  }

  if (use_beamspot_corr()) retrieve_beamspot_ = true;

  // LAr currents are required to calculate LAr lumi
  if (retrieve_lumi_LAr_) retrieve_currents_ = true;

  TRY( ReadChannels() )

  if (retrieve_lumi_LAr_) {
    TRY( ReadCalibrations(&channel_calibrations_, primary_calibrations_filepath()) )
    if (use_baseline_subtraction_from_fit()) {
      for (auto& channel: channel_calibrations_) {
        channel.second.intercept = 0.0;
      }
    }
  }

  return Void();
}

// Reads in calibration values for each of the channels being used (those
//   read in with ReadChannels)
Expected<Void> Analysis::ReadCalibrations(
    map<string, Analysis::ChannelCalibration>* channel_calibrations,
    string primary_calibrations_filepath)
{
  auto this_func_name = "Analysis::ReadCalibrations";
  if (channel_calibrations == nullptr) {
    return make_unexpected(make_shared<Error::Nullptr>(this_func_name));
  }
  for (auto &channel: *channel_calibrations) {
    std::ifstream calibrations_file(primary_calibrations_filepath);
    if (!calibrations_file) {
      return make_unexpected(make_shared<Error::File>(primary_calibrations_filepath, this_func_name));
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

  std::ifstream channels_list_file(channels_list_filepath());
  if (!channels_list_file) {
    return make_unexpected(make_shared<Error::File>(channels_list_filepath(), this_func_name));
  }

  string channel_name;
  while ( getline(channels_list_file, channel_name) ) {
    if (channel_name[0] == '#') {
      if (verbose()) cout << "Skipping channel " << channel_name << endl;
      continue;
    }
    channel_calibrations_.emplace_hint(
        end(channel_calibrations_),
        channel_name,
        ChannelCalibration{0.0, 0.0}
    );
  }
  channels_list_file.close();

  return Void();
}

// Control flow for the analysis of samples and creation of plots
Expected<Void> Analysis::RunAnalysis()
{
  auto run_names_vec = VectorFromFile(run_list_dir());
  RETURN_IF_ERR( run_names_vec )

  // Iterate over samples
  map<string, Run> runs_data;

  for (const auto &run_name: *run_names_vec) {
    auto this_run = Run(run_name, this);
    TRY_CONTINUE( this_run.Init() )
    TRY_CONTINUE( this_run.CreateSingleRunPlots() )

    if (retrieve_lumi_LAr_) {
      TRY_CONTINUE( this_run.CalcLArLumi() )
    }
    if (retrieve_mu_LAr_) {
      TRY_CONTINUE( this_run.CalcLArMu() )
    }
    if (do_benedetto()) {
      TRY_CONTINUE( this_run.CreateBenedettoOutput() )
    }

    runs_data.emplace_hint(end(runs_data),  //change this to cend in C++14
                           string(run_name), std::move(this_run));
    // DELET THIS
    auto node = "mode_options.lumi_current.";
    LumiCurrentPlotOptions plot_options(params_, node);
    plot_options.run_name(string("multi_to_")+run_name);
    for (const auto& channel: channel_calibrations()) {
      auto this_channel_name = channel.first;
      auto lumi_current_points = GenerateLumiVsCurrentPoints(runs_data, this_channel_name);
      auto is_verbose = verbose();
      lumi_current_points.erase(
          std::remove_if(
              lumi_current_points.begin(),
              lumi_current_points.end(),
              [is_verbose] (auto point) {
                auto lumi = point[0];
                auto current = point[1];
                if (is_verbose && current < 0) {
                  // Placeholder current representing an absent value is -9999.0
                  if (fabs(current + 9999.0) < 1.0) {
                    cout << "No current for this LB" << endl;
                  }
                  else {
                    cout << "Negative current. Current == " << current << ", lumi == " << lumi << endl;
                  }
                }
                return lumi < gOflLumiCutoff || current < gLArCurrentCutoff;
              }
              ),
          lumi_current_points.end());
      plot_options.channel_name(std::move(this_channel_name));
      plot_options.title("Runs "+runs_data.begin()->first+" - "+runs_data.rbegin()->first+", Channel "+this_channel_name);
      auto this_channel_fit_results = Plotter::PlotLumiCurrent(lumi_current_points, plot_options);
      CONTINUE_IF_ERR(this_channel_fit_results)
    }
  }

  CreateAllRunPlots(runs_data);

  return Void();
}

const map<string, int>& Analysis::n_bunches_from_file()
{ 
    if (!n_bunches_from_file_) {
      string n_bunches_file_path_node = "input_filepaths.n_bunches";
      auto n_bunches_file_path = params_.get<string>(n_bunches_file_path_node);
      if (!n_bunches_file_path || n_bunches_file_path->empty()) {
        n_bunches_from_file_ = map<string, int>();
      }
      else {
        auto n_bunches_file = JSONReader(*n_bunches_file_path);
        n_bunches_from_file_ = n_bunches_file.get<map<string, int>, int>("n_bunches");
        if (!n_bunches_from_file_) n_bunches_from_file_ = map<string, int>();
      }
    }
    return *n_bunches_from_file_;
}

const map<string, vector<int>>& Analysis::custom_LB_bounds()
{ 
    if (!custom_LB_bounds_) {
      custom_LB_bounds_ = map<string, vector<int>>();
      auto custom_LB_bounds_file_path_node = "input_filepaths.custom_LB_bounds";
      auto custom_LB_bounds_file_path = params_.get<string>(custom_LB_bounds_file_path_node);
      if (custom_LB_bounds_file_path && !custom_LB_bounds_file_path->empty()) {
        auto custom_LB_bounds_file = JSONReader(*custom_LB_bounds_file_path);
        auto custom_LB_bounds_from_file = custom_LB_bounds_file.get<map<string, vector<int>>, vector<int>, int>("custom_LB_bounds");
        if (custom_LB_bounds_from_file) {
          auto verified = VerifyLBBounds(*custom_LB_bounds_from_file);
          if (verified.valid()) {
            custom_LB_bounds_ = *custom_LB_bounds_from_file;
          }
          else {
            cout << "Error: one or more LB bounds settings are illegal" << endl;
          }
        }
      }
    }
    return *custom_LB_bounds_;
}

Detector::Name Analysis::detector() {
  if (!detector_) {
    detector_ = Detector::Name::FCal;
    auto detector_node = "detector";
    auto detector_str = params_.get<string>(detector_node);
    if (detector_str) {
      auto detector_name = Detector::FromString(std::move(*detector_str));
      if (detector_name.valid()) {
        detector_ = *detector_name;
      }
      else {
        cout << "Invalid detector name from param tree. Defaulting to FCal." << endl;
      }
    }
  }
  return *detector_;
}

#define LAZY_LOAD_AND_RETURN(member, type, node_name) \
  if (!member) { \
    member = params_.get<type>(node_name); \
  } \
  return *member;

bool Analysis::verbose() {
  if (!verbose_) {
    verbose_ = params_.verbose(); 
  }
  return *verbose_;
}

double Analysis::f_rev() {
  LAZY_LOAD_AND_RETURN(f_rev_, double, "lumi_calculation.f_rev")
}
double Analysis::x_sec() {
  LAZY_LOAD_AND_RETURN(x_sec_, double, "lumi_calculation.x_sec")
}

int Analysis::ref_run_number() {
  LAZY_LOAD_AND_RETURN(ref_run_number_, int, "ref_run_number")
}
double Analysis::anchoring_factor_A() {
  LAZY_LOAD_AND_RETURN(anchoring_factor_A_, double, "anchoring_factors."+std::to_string(ref_run_number())+".A")
}
double Analysis::anchoring_factor_C() {
  LAZY_LOAD_AND_RETURN(anchoring_factor_C_, double, "anchoring_factors."+std::to_string(ref_run_number())+".C")
}
double Analysis::anchoring_factor_Avg() {
  LAZY_LOAD_AND_RETURN(anchoring_factor_Avg_, double, "anchoring_factors."+std::to_string(ref_run_number())+".Avg")
}

string INPUT_FILEPATHS_NODE = "input_filepaths.";

string Analysis::primary_calibrations_filepath() {
  LAZY_LOAD_AND_RETURN(primary_calibrations_filepath_, string, INPUT_FILEPATHS_NODE+"primary_calibrations")
}
string Analysis::calibrations_dir() {
  LAZY_LOAD_AND_RETURN(calibrations_dir_, string, INPUT_FILEPATHS_NODE+"calibrations_dir")
}
string Analysis::channels_list_filepath() {
  LAZY_LOAD_AND_RETURN(channels_list_filepath_, string, INPUT_FILEPATHS_NODE+"channels_list")
}
string Analysis::pedestals_dir() {
  LAZY_LOAD_AND_RETURN(pedestals_dir_, string, INPUT_FILEPATHS_NODE+"pedestals")
}
string Analysis::trees_dir() {
  LAZY_LOAD_AND_RETURN(trees_dir_, string, INPUT_FILEPATHS_NODE+"trees")
}
string Analysis::currents_dir() {
  LAZY_LOAD_AND_RETURN(currents_dir_, string, INPUT_FILEPATHS_NODE+"currents")
}
string Analysis::run_list_dir() {
  LAZY_LOAD_AND_RETURN(run_list_dir_, string, INPUT_FILEPATHS_NODE+"run_list")
}

// Lazy load macro doesn't work here because get is taking two type parameters
vector<string> Analysis::modes() {
  if (!modes_) {
    modes_ = params_.get<vector<string>, string>("modes");
  }
  return *modes_;
}

string Analysis::reference_lumi_algo() {
  LAZY_LOAD_AND_RETURN(reference_lumi_algo_, string, "reference_lumi_algo")
}

bool Analysis::apply_LUCID_mu_corr() {
  LAZY_LOAD_AND_RETURN(apply_LUCID_mu_corr_, bool, "apply_LUCID_mu_corr")
}

bool Analysis::use_beamspot_corr() {
  LAZY_LOAD_AND_RETURN(use_beamspot_corr_, bool, "beamspot_correction.use")
}

// Lazy load macro doesn't work here because get is taking two type parameters
vector<Float_t> Analysis::beamspot_corr_params() {
  if (!beamspot_corr_params_) {
    beamspot_corr_params_ = params_.get<vector<Float_t>, Float_t>("beamspot_correction.params");
  }
  return *beamspot_corr_params_;
}

bool Analysis::use_start_of_fill_pedestals() {
  LAZY_LOAD_AND_RETURN(use_start_of_fill_pedestals_, bool, "use_start_of_fill_pedestals")
}
bool Analysis::use_baseline_subtraction_from_fit() {
  LAZY_LOAD_AND_RETURN(use_baseline_subtraction_from_fit_, bool, "use_baseline_subtraction_from_fit")
}
bool Analysis::do_benedetto() {
  if (!do_benedetto_) {
    const auto& modes_vec = modes();
    do_benedetto_ = ( std::find(modes_vec.begin(), modes_vec.end(), "benedetto") != modes_vec.end() );
  }
  return *do_benedetto_;
}

string Analysis::benedetto_output_dir() {
  if (!benedetto_output_dir_) {
    benedetto_output_dir_ = *params_.get<string>("base_output_dir") +
      *params_.get<string>("mode_options.benedetto.output_dir");
  }
  return *benedetto_output_dir_;
}
