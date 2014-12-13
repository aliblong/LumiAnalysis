#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <sstream>
#include <tuple>
//#include <unordered_map>
#include <vector>

#include "TAxis.h"
#include "TBranch.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TGraph.h"
#include "THStack.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TProfile.h"
#include "TTree.h"

#include "cutoffs.h"
#include "fcal_region_data.h"
#include "fit_results.h"
#include "plotter.h"
#include "lumi_current_plot_options.h"
#include "mu_stab_plot_options.h"
#include "make_unique.h"

using std::cout;
using std::cerr;
using std::endl;

using std::string;
using std::vector;

typedef std::map<string, FitResults> FitResultsMap;

namespace {

enum class Sign {Positive, Negative};
enum class Axis {X, Y};

typedef std::tuple<FCalRegion, Axis, Sign> GeometricRegion;

struct GeometricRegionOrderer {
  bool operator() (const GeometricRegion &r1, const GeometricRegion &r2) const {
    auto region1 = std::get<0>(r1);
    auto region2 = std::get<0>(r2);
    auto axis1 = std::get<1>(r1);
    auto axis2 = std::get<1>(r2);
    auto sign1 = std::get<2>(r1);
    auto sign2 = std::get<2>(r2);

    if (region1 == region2) {
      if (axis1 == axis2) {
        if ( (sign1 == Sign::Positive) && (sign2 == Sign::Negative) ) {
          return true;
        } else {
          return false;
        }
      } else if (axis1 == Axis::X) {
        return true;
      } else {
        return false;
      }
    } else if (region1 == FCalRegion::A) {
      return true;
    } else if ( (region1 == FCalRegion::C) && (region2 == FCalRegion::Avg) ) {
      return true;
    } else {
      return false;
    }
  }
};

typedef std::map<GeometricRegion, Float_t, GeometricRegionOrderer>
          GeometricRegionMap;

string RegionToString(FCalRegion region) {
  switch (region) {
    case FCalRegion::A:
      return "A";
    case FCalRegion::C:
      return "C";
    case FCalRegion::Avg:
      return "Avg";
    default:
      cerr << "ERROR: in RegionToString(FCalRegion region)" << endl;
      return "ERROR";
  }
}

string AxisToString(Axis axis) {
  switch (axis) {
    case Axis::X:
      return "X";
    case Axis::Y:
      return "Y";
    default:
      cerr << "ERROR: in AxisToString(Axis axis)" << endl;
      return "ERROR";
  }
}

string SignToString(Sign sign) {
  switch (sign) {
    case Sign::Positive:
      return "+";
    case Sign::Negative:
      return "-";
    default:
      cerr << "ERROR: in SignToString(Sign sign)" << endl;
      return "ERROR";
  }
}

void DumpGeoRegion(const GeometricRegion &region) {
  cout << RegionToString(std::get<0>(region))
       << AxisToString(std::get<1>(region))
       << SignToString(std::get<2>(region)) << endl;
}

int WriteGeometricAnalysisToText(
      const std::map<string, Float_t> phi_slices_slopes_sums,
      const GeometricRegionMap &regions_slopes_sums,
      string run_name,
      string output_dir) {

  int err = system( ("mkdir -p "+output_dir).c_str() );

  string out_filepath = output_dir + run_name + ".dat";
  std::ofstream out_file(out_filepath);
  if (!out_file.is_open()) {
    cerr << "ERROR: in WriteGeometricAnalysisToText(...):" << endl;
    cerr << "\tcould not open file \'" << out_filepath << "\'" << endl;
    return 1;
  }

  for (const auto& this_phi_slice_slopes_sum: phi_slices_slopes_sums) {
    const auto &phi_slice_name = this_phi_slice_slopes_sum.first;
    const auto &slopes_sum = this_phi_slice_slopes_sum.second;

    out_file << phi_slice_name << ' ' << slopes_sum << '\n';
  }

  out_file << '\n';

  vector<FCalRegion> regions = {FCalRegion::A, FCalRegion::C, FCalRegion::Avg};
  vector<Axis> axes = {Axis::X, Axis::Y};
  vector<Sign> signs = {Sign::Positive, Sign::Negative};
  for (auto region: regions) {
    for (auto axis: axes) {
      for (auto sign: signs) {
        auto this_geometric_region = std::make_tuple(region, axis, sign);
        out_file << RegionToString(region) << AxisToString(axis)
                 << SignToString(sign) << ' '
                 << regions_slopes_sums.at(this_geometric_region) <<'\n';
      }
    }
  }

  out_file.close();
  return 0;
}

GeometricRegionMap InitRegionsSlopeSumsMap() {
  GeometricRegionMap result;
  std::vector<FCalRegion> regions = {FCalRegion::A,
                                     FCalRegion::C,
                                     FCalRegion::Avg};
  std::vector<Axis> axes = {Axis::X, Axis::Y};
  std::vector<Sign> signs = {Sign::Positive, Sign::Negative};

  for (auto region: regions) {
    for (auto axis: axes) {
      for (auto sign: signs) {
        result.insert(std::make_pair(std::make_tuple(region, axis, sign), 0.0));
        //DumpGeoRegion(std::make_tuple(region, axis, sign));
      }
    }
  }
  return result;
}

string StringFromInt(int i, unsigned width) {
  std::stringstream result;
  result.width(2);
  result.fill('0');
  result << i;
  return result.str();
}

int IntFromChar(char c) {
  int result = c - '0';
  return result;
}

Sign GetSign(Axis axis, string phi_slice_name) {
  auto phi_slice_number = atoi(phi_slice_name.substr(3,2).c_str());
  switch(axis) {
    case Axis::X:
      if ((phi_slice_number < 4) || (phi_slice_number > 11)) {
        return Sign::Positive;
      } else {
        return Sign::Negative;
      }
    case Axis::Y:
      if (phi_slice_number < 8) {
        return Sign::Positive;
      } else {
        return Sign::Negative;
      }
    default:
      cerr << "ERROR: in GetSign(...) - with Sign enum" << endl;
      return Sign::Positive;
  }
}

std::map<string, Float_t> InitPhiSlicesSlopeSumsMap() {
  std::map<string, Float_t> result;
  vector<string> phi_slice_names;
  for (unsigned i = 0; i < 16; ++i) {
    phi_slice_names.emplace_back("A1." + StringFromInt(i, 2));
    phi_slice_names.emplace_back("C1." + StringFromInt(i, 2));
  }
  for (const auto &phi_slice: phi_slice_names) {
    result.insert(std::make_pair(phi_slice, 0.0));
  }
  return result;
}

string GetPhiSliceFromID(char region_ID, char module_ID, char channel_ID) {
  auto module_ID_int = IntFromChar(module_ID);
  if ( (module_ID_int < 0) || (module_ID_int > 9) ) {
    cerr << "ERROR: in IntFromChar(char c) - module_ID input falls outside [0,9]"
         << endl;
    return "ERROR";
  }
  auto channel_ID_int = IntFromChar(channel_ID);
  if ( (channel_ID_int < 0) || (channel_ID_int > 9) ) {
    cerr << "ERROR: in IntFromChar(char c) - channel_ID input falls outside [0,9]"
         << endl;
    return "ERROR";
  }
  unsigned phi_slice_ID_int = (module_ID_int * 2) + (channel_ID_int / 4);
  string region_ID_str(1,region_ID);
  string phi_slice = region_ID_str + "1." + StringFromInt(phi_slice_ID_int, 2);
  return phi_slice;
}

string GetPhiSliceFromChannel(string channel_name) {
// This is vaguely hack-y & only works with a very specific channel-naming scheme

  FCalRegion this_region = FCalRegion::A;
  if (channel_name.at(1) == '1') {
    this_region = FCalRegion::C;
  }

  char region_ID;
  char module_ID;
  char channel_ID;
  switch (this_region) {
    case FCalRegion::A:
      region_ID = 'A';
      module_ID = channel_name.at(2);
      channel_ID = channel_name.at(4);
      break;
    case FCalRegion::C:
      region_ID = 'C';
      module_ID = channel_name.at(3);
      channel_ID = channel_name.at(5);
      break;
  }
  return GetPhiSliceFromID(region_ID, module_ID, channel_ID);
}

void FormatTHStack(const MuStabPlotOptions &plot_options,
                   THStack &stack) {
  Double_t ratio_max;
  Double_t ratio_min;

  if(plot_options.y_auto_range()) {
    ratio_max = stack.GetMaximum() * 1.2;
    ratio_min = stack.GetMinimum() * 0.8;
  } else {
    ratio_max = plot_options.y_max();
    ratio_min = plot_options.y_min();
  }

  stack.GetXaxis()->SetTimeDisplay(1);
  stack.GetXaxis()->SetTimeFormat("%m/%d");
  stack.GetXaxis()->SetTitle(plot_options.x_title().c_str());
  stack.GetYaxis()->SetTitle(plot_options.y_title().c_str());
  stack.SetMinimum(ratio_min);
  stack.SetMaximum(ratio_max);
}

vector<FCalRegionData> InitRegionsData(const MuStabPlotOptions &plot_options) {
// Creates a vector of plotting parameters for the A- and C- side <mu>
//   stability profiles, as well as the average

  vector<FCalRegionData> regions;
  regions.emplace_back(FCalRegion::A,
                       "mu_stab_A",
                       "FCal A1",
                       plot_options.marker_color_A(),
                       plot_options.marker_size_A(),
                       plot_options.marker_style_A());
  regions.emplace_back(FCalRegion::C,
                       "mu_stab_C",
                       "FCal C1",
                       plot_options.marker_color_C(),
                       plot_options.marker_size_C(),
                       plot_options.marker_style_C());
  /*
  regions.emplace_back(FCalRegion::Avg,
                       "mu_stab_Avg",
                       "FCal Avg",
                       plot_options.marker_color_Avg(),
                       plot_options.marker_size_Avg(),
                       plot_options.marker_style_Avg());
                       */
  return regions;
}

int PopulateTProfile(FCalRegion region,
                     const std::map<string, SingleRunData> &runs_data,
                     TProfile *profile) {
// Fills the TProfile for a given region using run data
  if (profile == nullptr) {
    cerr << "ERROR: nullptr passed to Plotter::PopulateTProfile(..., TProfile *)"
         << endl;
    return 1;
  }

  Float_t last_run_diff = 0.0;
  Int_t last_bin = -999; // may as well make use of the fact that it's signed
  for (const auto &run: runs_data) {
    const auto &run_name = run.first;
    const auto &this_run_data = run.second;
    const auto &timestamp = this_run_data.timestamp();

    auto this_bin = profile->GetXaxis()->FindBin(timestamp);
    if (this_bin == last_bin) {
      cerr << "Warning: bin collision for run " << run_name << endl;
    } else {
      last_bin = this_bin;
    }
    //cout << "########################" << endl;
    //cout << "Run " << run.first << ", timestamp = " << timestamp << endl;
    //unsigned ratios_counter = 0;
    //Float_t ratios = 0.0;
    auto n_events = this_run_data.lumi_BCM_.size();

    for (unsigned iEvent = 0; iEvent < n_events; ++iEvent) {
      auto event_lumi_BCM = this_run_data.lumi_BCM_.at(iEvent);
      //cout << "BCM lumi: " << event_lumi_BCM << endl;
      // Skip events where BCM lumi ~= 0
      if (event_lumi_BCM < gBCMLumiCutoff) continue;

      Float_t event_lumi_FCal;
      switch (region) {
        case FCalRegion::A:
          event_lumi_FCal = this_run_data.lumi_FCal_A_.at(iEvent);
          break;
        case FCalRegion::C:
          event_lumi_FCal = this_run_data.lumi_FCal_C_.at(iEvent);
          break;
        case FCalRegion::Avg:
          event_lumi_FCal = ( this_run_data.lumi_FCal_A_.at(iEvent) +
                              this_run_data.lumi_FCal_C_.at(iEvent) ) / 2;
          // Skip events where A- OR C-side FCal lumi ~= 0
          if (this_run_data.lumi_FCal_A_.at(iEvent) < gFCalLumiCutoff ||
              this_run_data.lumi_FCal_C_.at(iEvent) < gFCalLumiCutoff) continue;
          break;
        default:
          cerr << "Error using FCalRegion enum" << endl;
          event_lumi_FCal = -999.0;
      }
      //cout << "FCal lumi: " << event_lumi_FCal << endl;
      if (event_lumi_FCal < gFCalLumiCutoff) continue;
      // Plot the percent difference between FCal lumi & BCM lumi
      Float_t this_event_diff = ((event_lumi_FCal/event_lumi_BCM) - 1)*100;
      if (this_event_diff < gPercentDiffMin ||
          this_event_diff > gPercentDiffMax) continue;
      /*
      if (run_name == "208642") {
        cout << event_lumi_FCal << endl;
        cout << event_lumi_BCM << endl;
        cout << this_event_diff << endl;
        cout << "#########################" << endl;
      }
      */
      profile->Fill(timestamp, this_event_diff);
      //ratios += event_lumi_BCM/event_lumi_FCal;
      //++ratios_counter;
      /*
      if (run_name == "202668" ||
          run_name == "202712" ||
          run_name == "208179" ||
          run_name == "208642" ||
          run_name == "211620" ||
          run_name == "211670") {
        TLine
      */
    }

    /*
    Float_t this_run_diff = profile->GetBinContent(profile->GetXaxis()->FindBin(timestamp));
    if ((abs(this_run_diff - last_run_diff) > 0.5) && (this_run_diff > 0)) {
      cout << "######################################" << endl;
      cout << "Anomalous run: " << run_name << endl;
      cout << "Bin value: " << this_run_diff << endl;
      cout << "Number of coll. bunches: " << this_run_data.nCollisions_ << endl;
    } else {
      last_run_diff = this_run_diff;
    }
    */
  }

  return 0;
}

void ScaleVector(vector<Float_t> &vec, Float_t factor) {
  for (auto &element: vec) {
    element *= factor;
  }
}

void SetAxisAutoRange(const vector<Float_t> &points, TAxis *axis) {
  auto max = *(std::max_element(points.begin(), points.end()));
  auto min_nonzero = *(std::min_element(points.begin(), points.end()));
  axis->SetRangeUser(0.8*min_nonzero, 1.2*max);
}

string TruncateFloatToString(Float_t value, unsigned n_decimal_places) {
  std::stringstream stream;
  stream << std::fixed << std::setprecision(n_decimal_places) << value;
  return stream.str();
}

void FormatLegend(const FitResults &fit_results, TPaveText &fit_legend) {
  Float_t goodness_of_fit = fit_results.chi_squared / fit_results.nDoF;
  string GoF = TruncateFloatToString(goodness_of_fit, 2);
  string slope = TruncateFloatToString(fit_results.slope, 3);
  string slope_err = TruncateFloatToString(fit_results.slope_err, 3);
  string intercept = TruncateFloatToString(fit_results.intercept, 2);
  string intercept_err = TruncateFloatToString(fit_results.intercept_err, 2);
  fit_legend.AddText( ("#chi^{2}/nDoF = "+GoF).c_str() );
  fit_legend.AddText( ("slope = "+slope+" #pm "+slope_err).c_str() );
  fit_legend.AddText( ("constant = "+intercept+" #pm "+intercept_err).c_str() );
  fit_legend.SetX1NDC(0.15);
  fit_legend.SetX2NDC(0.4);
  fit_legend.SetY1NDC(0.65);
  fit_legend.SetY2NDC(0.88);
  fit_legend.SetTextSize(0.035);
  fit_legend.SetOption("");
  fit_legend.SetBorderSize(0);
  fit_legend.SetFillColor(0);
}

/* Unused functions
template<typename T>
struct BinaryAverage {
  Float_t operator() (const T &T1, const T &T2) { return (T1 + T2) / 2.0; }
};

void DumpVector(const vector<Float_t> &vec) {
  for (auto element: vec) cout << element << endl;
}

bool MinNonZero(Float_t i, Float_t j) {
  Float_t min_allowable = 0.1;
  if (j < min_allowable) {
    return 1;
  } else if (i < min_allowable) {
    return 0;
  } else {
    return i < j;
  }
  return i < j;
}
*/

}

int Plotter::PlotLumiCurrent(const vector<Float_t> &lumi_arg,
                             const vector<Float_t> &current_arg,
                             string run_name,
                             string channel_name,
                             const LumiCurrentPlotOptions &plot_options,
                             string output_dir,
                             FitResults &fit_results) {
  TCanvas canvas;
  canvas.cd();

  vector<Float_t> lumi;
  vector<Float_t> current;
  auto num_points_arg = lumi_arg.size();
  Float_t epsilon = 0.00000000001;

  // Filter out events with lumi or current of 0
  for (unsigned i=0; i < num_points_arg; ++i) {
    if (lumi_arg.at(i) > epsilon && current_arg.at(i) > epsilon) { //&&
        //!(run_name == "203934" && i == num_points_arg - 1)) {
      lumi.push_back(lumi_arg.at(i));
      current.push_back(current_arg.at(i));
    }
  }

  // Don't bother running the scaling transform if scale_factor = 1.0
  if (abs(plot_options.x_scale() - 1) < epsilon) {
    ScaleVector(lumi, plot_options.x_scale());
  }
  if (abs(plot_options.y_scale() - 1) < epsilon) {
    ScaleVector(current, plot_options.y_scale());
  }

  // Use C-style dynamic arrays since ROOT doesn't support vectors
  auto num_points = lumi.size();
  Float_t *lumi_arr = new Float_t[num_points];
  Float_t *current_arr = new Float_t[num_points];
  for (unsigned iPoint=0; iPoint < num_points; ++iPoint) {
    lumi_arr[iPoint] = lumi.at(iPoint);
    current_arr[iPoint] = current.at(iPoint);
  }

  TGraph graph(num_points, lumi_arr, current_arr);

  // Plot formatting
  graph.SetMarkerColor(plot_options.marker_color());
  graph.SetMarkerStyle(plot_options.marker_style());
  graph.SetMarkerSize(plot_options.marker_size());

  string graph_title = "Run "+run_name+", Channel "+channel_name;
  graph.SetTitle( graph_title.c_str() );

  // Axes formatting
  // Must call Draw before doing anything with the axes
  graph.Draw(plot_options.draw_options().c_str());

  if (plot_options.x_auto_range()) {
    SetAxisAutoRange(lumi, graph.GetXaxis());
  } else {
    graph.GetXaxis()->SetRangeUser(plot_options.x_min(), plot_options.x_max());
  }
  if (plot_options.y_auto_range()) {
    SetAxisAutoRange(current, graph.GetYaxis());
  } else {
    graph.GetYaxis()->SetRangeUser(plot_options.y_min(), plot_options.y_max());
  }

  graph.GetXaxis()->SetTitle(plot_options.x_title().c_str());
  graph.GetYaxis()->SetTitle(plot_options.y_title().c_str());

  // Write fit results and format fit legend
  TPaveText fit_legend;
  if (plot_options.do_fit()) {
    graph.Fit("pol1", plot_options.fit_options().c_str());
    TF1 *fit = graph.GetFunction("pol1");
    fit->SetLineColor(plot_options.fit_line_color());
    fit->SetLineWidth(plot_options.fit_line_width());

    fit_results.FromFit(fit, plot_options);

    if (plot_options.fit_show_legend()) {
      FormatLegend(fit_results, fit_legend);
    }
  }

  //graph.Draw(plot_options.draw_options.c_str());
  fit_legend.Draw();

  canvas.Update();
  string write_dir = output_dir+"lumi_current/"+run_name+"/";
  int err_system = system( ("mkdir -p "+write_dir).c_str() );

  canvas.Print( (write_dir+channel_name+".png").c_str() );
  //TFile *this_file = TFile::Open("graph_file.root", "RECREATE");
  //graph.Write();

  //this_file->Close();
  delete [] lumi_arr;
  delete [] current_arr;

  return 0;
}

/*
int Plotter::GeometricAnalysisFromFitResultsTree(string run_name,
                                                 string output_dir) {
// Writes FCal current vs. BCM lumi fit parameters to a root file for a given
//   run. These parameters are used to derive the FCal lumi calibration.

  string read_dir = output_dir+"fit_results/";
  int err = system( ("mkdir -p "+read_dir).c_str() );

  TFile read_file((read_dir+run_name+".root").c_str(), "recreate");
  //string tree_name = "t_"+run_name;
  TTree *read_tree = static_cast<TTree*>("");

  string channel_name;
  read_tree->SetBranchAddress("channel_name", &channel_name);

  FitResults results_to_write;
  read_tree->SetBranchAddress("slope", &results_to_write.slope);
  read_tree->SetBranchAddress("slope_err", &results_to_write.slope_err);
  read_tree->SetBranchAddress("intercept", &results_to_write.intercept);
  read_tree->SetBranchAddress("intercept_err", &results_to_write.intercept_err);
  read_tree->SetBranchAddress("chi_squared", &results_to_write.chi_squared);
  read_tree->SetBranchAddress("nDoF", &results_to_write.nDoF);
  read_tree->SetBranchAddress("is_short", &results_to_write.is_short);
  read_tree->SetBranchAddress("calibration_slope",
                              &results_to_write.calibration_slope);
  read_tree->SetBranchAddress("calibration_intercept",
                              &results_to_write.calibration_intercept);

  for (const auto& this_channel_results: fit_results) {
    channel_name = this_channel_results.first;
    results_to_write = this_channel_results.second;

    tree.Fill();
  }

  tree.Write();

  return 0;
}
*/

int Plotter::GeometricAnalysisOfFitResults(
                     const FitResultsMap &fit_results,
                     string run_name,
                     string output_dir) {
// Writes FCal current vs. BCM lumi fit parameters to a root file for a given
//   run. These parameters are used to derive the FCal lumi calibration.

  auto phi_slices_slopes_sums = InitPhiSlicesSlopeSumsMap();
  for (const auto &this_channel_results: fit_results) {
    auto phi_slice = GetPhiSliceFromChannel(this_channel_results.first);
    phi_slices_slopes_sums.at(phi_slice) += this_channel_results.second.calibration_slope;
  }

  auto regions_slopes_sums = InitRegionsSlopeSumsMap();

  /*
  for (const auto &asdf: regions_slopes_sums) {
    DumpGeoRegion(asdf.first);
  }
  */

  vector<Axis> axes = {Axis::X, Axis::Y};
  vector<Sign> signs = {Sign::Positive, Sign::Negative};
  for (const auto &this_phi_slice_slopes_sum: phi_slices_slopes_sums) {
    const auto &phi_slice_name = this_phi_slice_slopes_sum.first;
    const auto &slopes_sum = this_phi_slice_slopes_sum.second;

    auto region = FCalRegion::A;
    if (phi_slice_name.at(0) == 'C') region = FCalRegion::C;

    for (auto axis: axes) {
      auto sign = GetSign(axis, phi_slice_name);
      GeometricRegion this_geometric_region = std::make_tuple(region, axis, sign);
      //DumpGeoRegion(this_geometric_region);
      regions_slopes_sums.at(this_geometric_region) += slopes_sum;
    }
  }

  for (auto axis: axes) {
    for (auto sign: signs) {
      auto region_A = std::make_tuple(FCalRegion::A, axis, sign);
      auto region_C = std::make_tuple(FCalRegion::C, axis, sign);
      auto region_Avg = std::make_tuple(FCalRegion::Avg, axis, sign);

      auto slopes_sum_A = regions_slopes_sums.at(region_A);
      auto slopes_sum_C = regions_slopes_sums.at(region_C);

      regions_slopes_sums.at(region_Avg) = (slopes_sum_A + slopes_sum_C) / 2;
    }
  }

  WriteGeometricAnalysisToText(phi_slices_slopes_sums,
                               regions_slopes_sums,
                               run_name,
                               output_dir);
  return 0;
}

int Plotter::WriteFitResultsToTree(const FitResultsMap &fit_results,
                                   string run_name,
                                   string output_dir) {
 // Writes FCal current vs. BCM lumi fit parameters to a root file for a given
 //   run. These parameters are used to derive the FCal lumi calibration.
 
   int err = system( ("mkdir -p "+output_dir).c_str() );
 
   TTree tree;
 
   string channel_name;
   tree.Branch("channel_name", &channel_name);
 
   FitResults results_to_write;
   tree.Branch("slope", &results_to_write.slope);
   tree.Branch("slope_err", &results_to_write.slope_err);
   tree.Branch("intercept", &results_to_write.intercept);
   tree.Branch("intercept_err", &results_to_write.intercept_err);
   tree.Branch("chi_squared", &results_to_write.chi_squared);
   tree.Branch("nDoF", &results_to_write.nDoF);
   tree.Branch("is_short", &results_to_write.is_short);
   tree.Branch("calibration_slope", &results_to_write.calibration_slope);
   tree.Branch("calibration_intercept", &results_to_write.calibration_intercept);
 
   for (const auto& this_channel_results: fit_results) {
     channel_name = this_channel_results.first;
     results_to_write = this_channel_results.second;
 
     tree.Fill();
   }

  TFile file((output_dir+run_name+".root").c_str(), "recreate");
  tree.Write();
  file.Close();

  return 0;
}

int Plotter::WriteCalibrationToText(const FitResultsMap &fit_results,
                                    string run_name,
                                    string output_dir) {
 // Writes the FCal lumi calibration data space-delimited in a plaintext file:
 //     
 //     ChannelName slope intercept

   int err = system( ("mkdir -p "+output_dir).c_str() );

   string out_filepath = output_dir + run_name + ".dat";
   std::ofstream out_file(out_filepath);
   if (!out_file.is_open()) {
     cerr << "ERROR: in WriteCalibrationToText(...):" << endl;
     cerr << "\tcould not open file \'" << out_filepath << "\'" << endl;
     return 1;
   }

   for (const auto& this_channel_results: fit_results) {
     const auto& result = this_channel_results.second;
 
     string channel_name = this_channel_results.first;
     Float_t slope = result.calibration_slope;
     Float_t intercept = result.calibration_intercept;
 
     out_file << channel_name << ' ' << slope << ' ' << intercept << '\n';
   }

   out_file.close();
   return 0;
 }

int Plotter::PlotMuStability(const std::map<string, SingleRunData> &runs_data,
                             const MuStabPlotOptions &plot_options,
                             string output_dir) {
// Plots TProfiles of <mu> time stability for select FCal regions (A-side,
//   C-side, average between them, etc.) all on the same canvas.

  TCanvas canvas;
  canvas.cd();
  THStack stack("stack","");

  TLegend legend(0.62,0.75,0.87,0.82);
  legend.SetFillColor(0);
  legend.SetBorderSize(1);
  legend.SetNColumns(2);

  // Plots and root files are saved in separate directories
  string write_dir = output_dir+"mu_stability/";
  string rootfiles_write_dir = write_dir+"histo_files/";
  int err_mkdir1 = system( ("mkdir -p "+write_dir).c_str() );
  int err_mkdir2 = system( ("mkdir -p "+rootfiles_write_dir).c_str() );

  // This vector holds the plotting parameters for each region
  auto regions = InitRegionsData(plot_options);

  // Binning is done with unix timestamp ranges
  UInt_t low_bin;
  UInt_t high_bin;
  UInt_t n_bins;
  if (plot_options.x_auto_range()) {
    low_bin = runs_data.begin()->second.timestamp();
    high_bin = runs_data.end()->second.timestamp();
    n_bins = 4320;
  } else {
    low_bin = plot_options.low_bin();
    high_bin = plot_options.high_bin();
    n_bins = plot_options.n_bins();
  }

  // Creates a TProfile for each region and adds it to the stack
  for (auto &region: regions) {

    // TProfiles must be added to the stack via owning pointers
    region.plot_ = std::make_unique<TProfile>(region.plot_name_.c_str(),
                                region.plot_title_.c_str(),
                                n_bins,
                                low_bin,
                                high_bin);

    int err_populate = PopulateTProfile(region.region_name_,
                                        runs_data,
                                        region.plot_.get());
    if (err_populate) {
      return 1;
    }

    // Profiles are saved individually in root files
    TFile file( (rootfiles_write_dir+region.plot_name_+".root").c_str(),
                "RECREATE" );
    region.plot_->Write();
    file.Close();

    region.plot_->SetMarkerColor(region.marker_color_);
    region.plot_->SetMarkerSize(region.marker_size_);
    region.plot_->SetMarkerStyle(region.marker_style_);

    stack.Add(region.plot_.get(), "AP");
    legend.AddEntry(region.plot_.get(), region.plot_title_.c_str(), "p");
  }

  // Draw must be called before Format! ROOT's fault, not mine
  stack.Draw(plot_options.draw_options().c_str());
  FormatTHStack(plot_options, stack);

  legend.Draw();

  std::vector<string> interest_samples = {"202668",
                                          "202712",
                                          "208179",
                                          "208642",
                                          "211620",
                                          "211670"};

  std::vector<TLine> interest_sample_markers;
  for (const auto &sample: interest_samples) {
    auto interest_run = runs_data.find(sample);
    if (interest_run != runs_data.end()) {
      auto timestamp = interest_run->second.timestamp();
      interest_sample_markers.emplace_back(timestamp, plot_options.y_min(),
                                           timestamp, plot_options.y_max());
    }
  }

  for (auto &sample_marker: interest_sample_markers) {
    sample_marker.SetLineColor(2);
    sample_marker.SetLineWidth(1);
    sample_marker.Draw();
  }

  TLatex label_ATLAS;
  label_ATLAS.SetNDC(true);
  label_ATLAS.SetTextSize(0.040);
  label_ATLAS.SetTextFont(72);
  label_ATLAS.DrawLatex(0.65,0.84,"ATLAS");
  label_ATLAS.SetTextFont(42);
  label_ATLAS.DrawLatex(0.65+0.10,0.84,"Internal");
  label_ATLAS.SetTextSize(0.030);

  canvas.Update();
  canvas.Print( (write_dir+"mu_stability.pdf").c_str() );

  return 0;
}
