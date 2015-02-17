#include "plotter.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <sstream>
#include <tuple>
#include <vector>

#include "TAxis.h"
#include "TBranch.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "THStack.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TProfile.h"
#include "TTree.h"

#include "boost/expected/expected.hpp"

#include "ATLAS_label_options.h"
#include "cutoffs.h"
#include "error.h"
#include "fcal_region_data.h"
#include "fcal_region.h"
#include "fit_results.h"
#include "lumi_current_plot_options.h"
#include "mu_stab_plot_options.h"
#include "point.h"
#include "util.h"
#include "void.h"

using std::cout;
using std::cerr;
using std::endl;

using std::string;
using std::vector;

using std::make_shared;

using boost::expected;
using boost::make_expected;
using boost::make_unexpected;

using FCalRegion::ZSide;
using FCalRegion::Axis;
using FCalRegion::Sign;

using Error::Expected;

typedef vector<Float_t> VectorF;
typedef std::map<string, FitResults> FitResultsMap;

namespace {

typedef std::map<FCalRegion::ModuleHalf, Float_t>
          ModuleHalfMap;

ModuleHalfMap InitRegionsSlopeSumsMap() {
  ModuleHalfMap result;
  auto module_half_set = FCalRegion::CreateModuleHalfSet();
  for (const auto& region: module_half_set) {
    result.insert(std::make_pair(region, 0.0));
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

Expected<Void> WriteGeometricAnalysisToText(
      const std::map<string, Float_t>& phi_slices_slopes_sums,
      const ModuleHalfMap& regions_slopes_sums,
      const LumiCurrentPlotOptions& options) {

  // For error reporting. Replace with reflection framework?
  const char *this_func_name = "WriteGeometricAnalysisToText";

  auto output_dir = options.geometric_results_dir();
  TRY( Util::mkdir(output_dir) )

  const string out_filepath = output_dir + options.run_name() + ".dat";
  std::ofstream out_file(out_filepath);
  if (!out_file.is_open()) {
    return make_unexpected(make_shared<Error::File>(out_filepath, this_func_name));
  }

  for (const auto& this_phi_slice_slopes_sum: phi_slices_slopes_sums) {
    const auto& phi_slice_name = this_phi_slice_slopes_sum.first;
    const auto& slopes_sum = this_phi_slice_slopes_sum.second;

    out_file << phi_slice_name << ' ' << slopes_sum << '\n';
  }

  out_file << '\n';

  auto regions = FCalRegion::CreateModuleHalfSet();
  for (const auto& region: regions) {
    auto z_side = std::get<0>(region);
    auto axis = std::get<1>(region);
    auto sign = std::get<2>(region);
    out_file << FCalRegion::ToString(z_side)
             << FCalRegion::ToString(axis)
             << FCalRegion::ToString(sign)
             << ' '
             << regions_slopes_sums.at(region)
             <<'\n';
  }

  out_file.close();
  return Void();
}

std::map<string, Float_t> InitPhiSlicesSlopeSumsMap() {
  std::map<string, Float_t> result;
  vector<string> phi_slice_names;
  for (unsigned i = 0; i < 16; ++i) {
    phi_slice_names.emplace_back("A1." + StringFromInt(i, 2));
    phi_slice_names.emplace_back("C1." + StringFromInt(i, 2));
  }
  for (const auto& phi_slice: phi_slice_names) {
    result.insert(std::make_pair(phi_slice, 0.0));
  }
  return result;
}

void FormatTHStack(const MuStabPlotOptions& plot_options,
                   THStack& stack) {
  Double_t ratio_max;
  Double_t ratio_min;

  if(plot_options.y_auto_range()) {
    ratio_max = stack.GetMaximum() * 1.2;
    ratio_min = stack.GetMinimum() * 0.8;
  }
  else {
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

vector<FCalRegionData> InitRegionsData(const MuStabPlotOptions& plot_options) {
// Creates a vector of plotting parameters for the A- and C- side <mu>
//   stability profiles, as well as the average

  vector<FCalRegionData> regions;
  regions.emplace_back(FCalRegion::ZSide::A,
                       "mu_stab_A",
                       "FCal A1",
                       plot_options.marker_color_A(),
                       plot_options.marker_size_A(),
                       plot_options.marker_style_A());
  regions.emplace_back(FCalRegion::ZSide::C,
                       "mu_stab_C",
                       "FCal C1",
                       plot_options.marker_color_C(),
                       plot_options.marker_size_C(),
                       plot_options.marker_style_C());
  /*
  regions.emplace_back(FCalRegion::ZSide::Avg,
                       "mu_stab_Avg",
                       "FCal Avg",
                       plot_options.marker_color_Avg(),
                       plot_options.marker_size_Avg(),
                       plot_options.marker_style_Avg());
                       */
  return regions;
}

// Fills the TProfile for a given side using run data
Expected<Void> PopulateTProfile(
    FCalRegion::ZSide side,
    const std::map<string, SingleRunData>& runs_data,
    TProfile *profile)
{
  auto func_name = "PopulateTProfile";
  if (profile == nullptr) {
    return make_unexpected(make_shared<Error::Nullptr>(func_name));
  }

  Int_t last_bin = -999;
  for (const auto& run: runs_data) {
    const auto& run_name = run.first;
    const auto& this_run_data = run.second;
    auto timestamp = this_run_data.timestamp();
    const auto& this_run_lumi_BCM = this_run_data.lumi_BCM();
    const auto& this_run_lumi_FCal_A = this_run_data.lumi_FCal_A();
    const auto& this_run_lumi_FCal_C = this_run_data.lumi_FCal_C();

    assert(this_run_lumi_BCM.size() == this_run_lumi_FCal_A.size());
    assert(this_run_lumi_FCal_A.size() == this_run_lumi_FCal_C.size());
    /*
    cout << "*******************************" << endl;
    cout << run_name << endl;
    cout << this_run_lumi_BCM.at(0) << endl;
    cout << this_run_lumi_FCal_A.at(0) << endl;
    cout << this_run_lumi_FCal_C.at(0) << endl;
    */

    auto this_bin = profile->GetXaxis()->FindBin(timestamp);
    if (this_bin == last_bin) {
      cerr << "Warning: bin collision for run " << run_name << endl;
    }
    else {
      last_bin = this_bin;
    }

    auto n_events = this_run_lumi_BCM.size();

    for (unsigned iEvent = 0; iEvent < n_events; ++iEvent) {
      auto event_lumi_BCM = this_run_lumi_BCM[iEvent];
      if (event_lumi_BCM < gBCMLumiCutoff) continue;

      Float_t event_lumi_FCal = 0.0;
      if (side == ZSide::A) {
        event_lumi_FCal = this_run_lumi_FCal_A[iEvent];
      }
      else if (side == ZSide::C) {
        event_lumi_FCal = this_run_lumi_FCal_C[iEvent];
      }
      else { // (side == ZSide::Both) {
        event_lumi_FCal = ( this_run_lumi_FCal_A[iEvent] +
                            this_run_lumi_FCal_C[iEvent] ) / 2;
        if (this_run_lumi_FCal_A[iEvent] < gFCalLumiCutoff ||
            this_run_lumi_FCal_C[iEvent] < gFCalLumiCutoff) continue;
      }
      if (event_lumi_FCal < gFCalLumiCutoff) continue;
      Float_t this_event_diff = ((event_lumi_FCal/event_lumi_BCM) - 1)*100;
      //if (this_event_diff < gPercentDiffMin ||
      //    this_event_diff > gPercentDiffMax) continue;

      profile->Fill(timestamp, this_event_diff);
    }

  }

  return Void();
}

void SetAxisAutoRange(TAxis *axis) {
  auto max = axis->GetXmax();
  auto min = axis->GetXmin();
  axis->SetRangeUser(0.8*min, 1.2*max);
}

string TruncateFloatToString(Float_t value, unsigned n_decimal_places) {
  std::stringstream stream;
  stream << std::fixed << std::setprecision(n_decimal_places) << value;
  return stream.str();
}

void FormatLegend(TPaveText& fit_legend, const FitResults& fit_results) {
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

void ApplyLCPlotOptionsToGraph(TGraphErrors& graph,
                               const LumiCurrentPlotOptions& plot_options) {
  graph.SetMarkerColor(plot_options.marker_color());
  graph.SetMarkerStyle(plot_options.marker_style());
  graph.SetMarkerSize(plot_options.marker_size());

  graph.GetXaxis()->SetRangeUser(plot_options.x_min(), plot_options.x_max());
  graph.GetYaxis()->SetRangeUser(plot_options.y_min(), plot_options.y_max());

  graph.GetXaxis()->SetTitle(plot_options.x_title().c_str());
  graph.GetYaxis()->SetTitle(plot_options.y_title().c_str());
}

TLatex DrawATLASLabel(const ATLASLabelOptions& options) {
  TLatex label;
  label.SetNDC(options.use_NDC());
  label.SetTextSize(options.text_size());
  label.SetTextFont(72);
  label.DrawLatex(options.x(),options.y(),"ATLAS");
  label.SetTextFont(42);
  label.DrawLatex(options.x()+options.subheading_offset(),
                  options.y(),
                  options.subheading().c_str());
  //label.SetTextSize(options.text_size());

  return label;
}

void FilterOutliers(TGraphErrors& graph, const TF1& fit,
                    const LumiCurrentPlotOptions& plot_options) {
  Double_t slope = fit.GetParameter(1);
  Double_t intercept = fit.GetParameter(0);

  Double_t x;
  Double_t y;

  Int_t nPoints = graph.GetN();

  vector<Double_t> to_be_removed;
  to_be_removed.reserve(nPoints);

  for (Int_t iPoint = 0; iPoint < nPoints; ++iPoint) {
    graph.GetPoint(iPoint, x, y);
    Double_t y_expected = slope*x + intercept;
    //cout << y_expected << ' ' << y << ' ' << graph.GetErrorY(iPoint) << endl;
    if ( std::abs(y - y_expected) > 20*graph.GetErrorY(iPoint)) {
      to_be_removed.push_back(iPoint);
    }
  }
  for (auto rit = to_be_removed.crbegin();
       rit != to_be_removed.crend();
       ++rit) {
    graph.RemovePoint(*rit);
    //cout << *rit << endl;
  }

  graph.Fit("f1", "BQS");

  //graph.Fit("pol1", plot_options.fit_options().c_str());
  //fit = graph.GetFunction("pol1");
}

void RecursiveFilterOutliers(TGraphErrors& graph, const TF1& fit,
                             const LumiCurrentPlotOptions& plot_options) {
  Int_t nPoints = graph.GetN();
  Int_t new_nPoints = 1;
  while (new_nPoints < nPoints && new_nPoints) {
    nPoints = graph.GetN();
    FilterOutliers(graph, fit, plot_options);
    new_nPoints = graph.GetN();
    //cout << nPoints << ' ' << new_nPoints << endl;
  }
}

void SetErrors(TGraphErrors& graph, Double_t x_rel_error, Double_t y_rel_error) {
  Int_t nPoints = graph.GetN();
  Double_t x;
  Double_t y;

  for (Int_t iPoint = 0; iPoint < nPoints; ++iPoint) {
    graph.GetPoint(iPoint, x, y);
    graph.SetPointError(iPoint, x_rel_error * x, y_rel_error * y);
  }
}

/* Unused functions
template<typename T>
struct BinaryAverage {
  Float_t operator() (const T& T1, const T& T2) { return (T1 + T2) / 2.0; }
};

void DumpVector(const VectorF& vec) {
  for (auto element: vec) cout << element << endl;
}

bool MinNonZero(Float_t i, Float_t j) {
  Float_t min_allowable = 0.1;
  if (j < min_allowable) {
    return 1;
  }
  else if (i < min_allowable) {
    return 0;
  }
  else {
    return i < j;
  }
  return i < j;
}
*/

}

Expected<FitResults> Plotter::PlotLumiCurrent(
    const VectorP<Float_t>& points,
    const LumiCurrentPlotOptions& options)
{
  const char *this_func_name = "PlotLumiCurrent";

  TCanvas canvas;
  canvas.cd();

  auto points_filtered = points;
  for (auto& point: points_filtered) {
    point[0] *= options.x_scale();
    point[1] *= options.y_scale();
  }

  points_filtered.erase(
      std::remove_if(
          points_filtered.begin(),
          points_filtered.end(),
          [] (auto point) {
            return point[0] < gFCalLumiCutoff || point[1] < gFCalCurrentCutoff;
          }
          ),
      points_filtered.end());

  // Use C-style dynamic arrays since ROOT doesn't support vectors
  auto num_points = points_filtered.size();
  Float_t *lumi_arr = new Float_t[num_points];
  Float_t *current_arr = new Float_t[num_points];
  for (unsigned iPoint=0; iPoint < num_points; ++iPoint) {
    lumi_arr[iPoint] = points_filtered[iPoint][0];
    current_arr[iPoint] = points_filtered[iPoint][1];
  }

  TGraphErrors graph(num_points, lumi_arr, current_arr);
  delete [] lumi_arr;
  delete [] current_arr;

  // Axes formatting
  // Must call Draw before doing anything with the axes
  graph.Draw(options.draw_options().c_str());
  SetErrors(graph, options.x_rel_error(), options.y_rel_error());

  ApplyLCPlotOptionsToGraph(graph, options);

  auto channel_name = options.channel_name();
  string graph_title = "Run "+options.run_name() +
                         ", Channel "+channel_name;
  graph.SetTitle( graph_title.c_str() );

  auto is_valid_channel = FCalRegion::IsValidChannel(channel_name);
  if (options.x_auto_range() || !is_valid_channel) SetAxisAutoRange(graph.GetXaxis());
  if (options.y_auto_range() || !is_valid_channel) SetAxisAutoRange(graph.GetYaxis());

  // Write fit results and format fit legend
  TPaveText fit_legend;
  FitResults fit_results;

  if (options.do_fit()) {
    TF1 fit("f1", "[0]+[1]*x", graph.GetXaxis()->GetXmin(), graph.GetXaxis()->GetXmax());

    if (options.fit_fix_intercept()) fit.FixParameter(0, 0.0);
    fit.SetLineColor(options.fit_line_color());
    fit.SetLineWidth(options.fit_line_width());

    graph.Fit("f1", (options.fit_options()+"BQS").c_str());

    //RecursiveFilterOutliers(graph, fit, options);

    fit_results.FromFit(fit, options);

    if (options.fit_show_legend()) {
      FormatLegend(fit_legend, fit_results);
    }
  }

  fit_legend.Draw();

  DrawATLASLabel(ATLASLabelOptions{"Internal"});

  canvas.Update();
  auto write_dir = options.plots_dir()+options.run_name()+"/";
  TRY( Util::mkdir(write_dir) )

  canvas.Print( (write_dir+options.channel_name()+".png").c_str() );
  return fit_results;
}

// Writes FCal current vs. BCM lumi fit parameters to a root file for a given
//   run. These parameters are used to derive the FCal lumi calibration.
Expected<Void> Plotter::GeometricAnalysisOfFitResults(
                     const FitResultsMap& fit_results,
                     const LumiCurrentPlotOptions& options)
{
  auto phi_slices_slopes_sums = InitPhiSlicesSlopeSumsMap();
  for (const auto& this_channel_results: fit_results) {
    auto phi_slice =
      FCalRegion::PhiSliceFromChannel(this_channel_results.first);
    RETURN_IF_ERR(phi_slice)
    phi_slices_slopes_sums.at(*phi_slice) += this_channel_results.second.slope;
  }

  auto regions_slopes_sums = InitRegionsSlopeSumsMap();

  for (const auto& this_phi_slice_slopes_sum: phi_slices_slopes_sums) {
    const auto& phi_slice_name = this_phi_slice_slopes_sum.first;
    const auto& slopes_sum = this_phi_slice_slopes_sum.second;

    auto region = FCalRegion::ZSide::A;
    if (phi_slice_name.at(0) == 'C') region = FCalRegion::ZSide::C;

    for (auto axis: FCalRegion::AXES) {
      auto sign = FCalRegion::SignFromAxisAndPhiSlice(axis, phi_slice_name);
      FCalRegion::ModuleHalf this_geometric_region = std::make_tuple(region, axis, sign);
      regions_slopes_sums.at(this_geometric_region) += slopes_sum;
    }
  }

  for (auto axis: FCalRegion::AXES) {
    for (auto sign: FCalRegion::SIGNS) {
      auto region_A = std::make_tuple(FCalRegion::ZSide::A, axis, sign);
      auto region_C = std::make_tuple(FCalRegion::ZSide::C, axis, sign);
      auto region_Avg = std::make_tuple(FCalRegion::ZSide::Both, axis, sign);

      auto slopes_sum_A = regions_slopes_sums.at(region_A);
      auto slopes_sum_C = regions_slopes_sums.at(region_C);

      regions_slopes_sums.at(region_Avg) = (slopes_sum_A + slopes_sum_C) / 2;
    }
  }

  return WriteGeometricAnalysisToText(phi_slices_slopes_sums,
                                      regions_slopes_sums,
                                      options);
}

// Writes FCal current vs. BCM lumi fit parameters to a root file for a given
//   run. These parameters are used to derive the FCal lumi calibration.
Expected<Void> Plotter::WriteFitResultsToTree(
    const FitResultsMap& fit_results,
    const LumiCurrentPlotOptions& options)
{
  const char *this_func_name = "WriteFitResultsToTree";

  auto output_dir = options.raw_fit_results_dir();
  TRY( Util::mkdir(output_dir) )
  TFile file((output_dir+options.run_name()+".root").c_str(),
              "recreate");
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

  tree.Write();
  file.Close();

  return Void();
}

// Writes the FCal lumi calibration data space-delimited in a plaintext file:
//
//     ChannelName slope intercept
Expected<Void> Plotter::WriteCalibrationToText(
    const FitResultsMap &fit_results,
    const LumiCurrentPlotOptions& options)
{
  const char *this_func_name = "WriteCalibrationToText";

  auto output_dir = options.calibration_results_dir();
  TRY( Util::mkdir(output_dir) )

  auto out_filepath = output_dir + options.run_name() + ".dat";
  std::ofstream out_file(out_filepath);
  if (!out_file.is_open()) {
    return make_unexpected(make_shared<Error::File>(out_filepath, this_func_name));
  }

  for (const auto& this_channel_results: fit_results) {
    const auto& result = this_channel_results.second;

    string channel_name = this_channel_results.first;
    Float_t slope = result.calibration_slope;
    Float_t intercept = result.calibration_intercept;

    out_file << channel_name << ' ' << slope << ' ' << intercept << '\n';
  }

  out_file.close();
  return Void();
}

// Plots TProfiles of <mu> time stability for select FCal regions (A-side,
//   C-side, average between them, etc.) all on the same canvas.
Expected<Void> Plotter::PlotMuStability(
    const std::map<string, SingleRunData> &runs_data,
    const MuStabPlotOptions &options)
{
  TCanvas canvas;
  canvas.cd();
  THStack stack("stack","");

  TLegend legend(0.62,0.75,0.87,0.82);
  legend.SetFillColor(0);
  legend.SetBorderSize(1);
  legend.SetNColumns(2);

  TRY( Util::mkdir(options.rootfiles_output_dir()) )

  // This vector holds the plotting parameters for each region
  auto regions = InitRegionsData(options);

  // Binning is done with unix timestamp ranges
  UInt_t low_bin;
  UInt_t high_bin;
  UInt_t n_bins;
  if (options.x_auto_range()) {
    low_bin = runs_data.begin()->second.timestamp();
    high_bin = runs_data.end()->second.timestamp();
    n_bins = 1000;
  }
  else {
    low_bin = options.low_bin();
    high_bin = options.high_bin();
    n_bins = options.n_bins();
  }

  // Creates a TProfile for each region and adds it to the stack
  for (auto &region: regions) {

    // TProfiles must be added to the THStack via owning pointers
    region.plot_ = std::make_unique<TProfile>(
        region.plot_name_.c_str(),
        region.plot_title_.c_str(),
        n_bins,
        low_bin,
        high_bin);

    TRY( PopulateTProfile(region.region_name_,
                          runs_data,
                          region.plot_.get()) )


    // Profiles are saved individually in root files
    TFile file( (options.rootfiles_output_dir()+region.plot_name_+".root").c_str(),
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
  stack.Draw(options.draw_options().c_str());
  FormatTHStack(options, stack);

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
      interest_sample_markers.emplace_back(timestamp, stack.GetYaxis()->GetXmin(),
                                           timestamp, stack.GetYaxis()->GetXmin());
    }
  }

  for (auto &sample_marker: interest_sample_markers) {
    sample_marker.SetLineColor(2);
    sample_marker.SetLineWidth(1);
    sample_marker.Draw();
  }

  DrawATLASLabel(ATLASLabelOptions{"Internal"});

  canvas.Update();
  canvas.Print( (options.base_output_dir()+"mu_stability.pdf").c_str() );

  return Void();
}

/*
int Plotter::PlotLumiTotalCurrent(const VectorF &lumi_arg,
                             const VectorF &current_arg_A,
                             const VectorF &current_arg_C,
                             string run_name,
                             const LumiCurrentPlotOptions &plot_options,
                             string output_dir) {
// This is a total hack since I'm short on time
  TCanvas canvas("c_lumi_tot_current", "placeholder", 1024, 768);
  canvas.cd();

  VectorF lumi;
  VectorF current_A;
  VectorF current_C;
  auto num_points_arg = lumi_arg.size();
  Float_t epsilon = 0.00000000001;

  // Filter out events with lumi or current of 0
  for (unsigned i=0; i < num_points_arg; ++i) {
    if (lumi_arg.at(i) > epsilon &&
        current_arg_A.at(i) > epsilon &&
        current_arg_C.at(i) > epsilon ) {
      lumi.push_back(lumi_arg.at(i));
      current_A.push_back(current_arg_A.at(i));
      current_C.push_back(current_arg_C.at(i));
    }
  }

  // Don't bother running the scaling transform if scale_factor = 1.0
  if (abs(plot_options.x_scale() - 1) < epsilon) {
    ScaleVector(lumi, plot_options.x_scale());
  }
  if (abs(plot_options.y_scale() - 1) < epsilon) {
    ScaleVector(current_A, plot_options.y_scale());
    ScaleVector(current_C, plot_options.y_scale());
  }

  // Use C-style dynamic arrays since ROOT doesn't support vectors
  auto num_points = lumi.size();
  Float_t *lumi_arr = new Float_t[num_points];
  Float_t *current_arr_A = new Float_t[num_points];
  Float_t *current_arr_C = new Float_t[num_points];
  for (unsigned iPoint=0; iPoint < num_points; ++iPoint) {
    lumi_arr[iPoint] = lumi.at(iPoint);
    current_arr_A[iPoint] = current_A.at(iPoint);
    current_arr_C[iPoint] = current_C.at(iPoint);
  }

  TGraphErrors graph_A(num_points, lumi_arr, current_arr_A);
  TGraphErrors graph_C(num_points, lumi_arr, current_arr_C);

  // Axes formatting
  // Must call Draw before doing anything with the axes
  graph_A.Draw(plot_options.draw_options().c_str());
  graph_C.Draw("PX");
  SetErrors(graph_A, plot_options.x_rel_error(), plot_options.y_rel_error());
  SetErrors(graph_C, plot_options.x_rel_error(), plot_options.y_rel_error());

  // Plot formatting
  ApplyLCPlotOptionsToGraph(graph_A, plot_options);
  ApplyLCPlotOptionsToGraph(graph_C, plot_options);
  graph_C.SetMarkerColor(3);
  //ApplyLCPlotOptionsToGraph(graph_C, plot_options);

  string graph_title = "Run "+run_name+", Sums";
  graph_A.SetTitle( graph_title.c_str() );

  if (plot_options.x_auto_range()) SetAxisAutoRange(graph_A.GetXaxis(), lumi);
  if (plot_options.y_auto_range()) SetAxisAutoRange(graph_A.GetYaxis(), current_A);
                                                        //  v~~~ not a typo!
  if (plot_options.y_auto_range()) SetAxisAutoRange(graph_C.GetYaxis(), current_A);

  TLatex label_ATLAS;
  DrawATLASLabel(label_ATLAS);

  canvas.Update();
  string write_dir = output_dir+"lumi_current/"+run_name+"/";
  int err_system = system( ("mkdir -p "+write_dir).c_str() );

  canvas.Print( (write_dir+"sums.png").c_str() );
  //TFile *this_file = TFile::Open("graph_A_file.root", "RECREATE");
  //graph.Write();

  //this_file->Close();
  delete [] lumi_arr;
  delete [] current_arr_A;
  delete [] current_arr_C;

  return 0;
}

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
