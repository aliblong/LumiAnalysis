#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <sstream>
#include <vector>

#include "TAxis.h"
#include "TBranch.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TGraph.h"
#include "THStack.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TProfile.h"
#include "TTree.h"

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

namespace {

vector<FCalRegionData> InitRegionsData(const MuStabPlotOptions &plot_options) {
// Creates a vector of plotting parameters for the A- and C- side <mu>
//   stability profiles, as well as the average

  vector<FCalRegionData> regions;
  regions.emplace_back(FCalRegion::A,
                       "mu_stab_A",
                       "FCal A",
                       plot_options.marker_color_A,
                       plot_options.marker_size_A,
                       plot_options.marker_style_A);
  regions.emplace_back(FCalRegion::C,
                       "mu_stab_C",
                       "FCal C",
                       plot_options.marker_color_C,
                       plot_options.marker_size_C,
                       plot_options.marker_style_C);
  regions.emplace_back(FCalRegion::Avg,
                       "mu_stab_Avg",
                       "FCal Avg",
                       plot_options.marker_color_Avg,
                       plot_options.marker_size_Avg,
                       plot_options.marker_style_Avg);
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

  for (const auto &run: runs_data) {
    const auto &this_run_data = run.second;
    const auto &timestamp = this_run_data.timestamp();

    unsigned ratios_counter = 0;
    float ratios = 0.0;
    auto n_events = this_run_data.lumi_BCM_.size();

    for (unsigned iEvent = 0; iEvent < n_events; ++iEvent) {
      auto event_lumi_BCM = this_run_data.lumi_BCM_.at(iEvent);
      // Skip events where BCM lumi ~= 0
      if (event_lumi_BCM < 0.1) continue;

      float event_lumi_FCal;
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
          if (this_run_data.lumi_FCal_A_.at(iEvent) < 0.1 ||
              this_run_data.lumi_FCal_C_.at(iEvent) < 0.1) continue;
          break;
        default:
          cerr << "Error using FCalRegion enum" << endl;
          event_lumi_FCal = -999.0;
      }
      if (event_lumi_FCal < 0.1) continue;
      // Plot the percent difference between FCal lumi & BCM lumi
      profile->Fill(timestamp, ((event_lumi_FCal/event_lumi_BCM) - 1)*100);
      ratios += event_lumi_BCM/event_lumi_FCal;
      ++ratios_counter;
    }
  }

  return 0;
}

void FormatTHStack(const MuStabPlotOptions &plot_options,
                   THStack &stack) {
  Double_t ratio_max;
  Double_t ratio_min;

  if(plot_options.y_auto_range) {
    ratio_max = stack.GetMaximum() * 1.2;
    ratio_min = stack.GetMinimum() * 0.8;
  } else {
    ratio_max = plot_options.y_max;
    ratio_min = plot_options.y_min;
  }

  stack.GetXaxis()->SetTimeDisplay(1);
  stack.GetXaxis()->SetTimeFormat("%m/%d");
  stack.GetXaxis()->SetTitle(plot_options.x_title.c_str());
  stack.GetYaxis()->SetTitle(plot_options.y_title.c_str());
  stack.SetMinimum(ratio_min);
  stack.SetMaximum(ratio_max);
}

template<typename T>
struct BinaryAverage {
  float operator() (const T &T1, const T &T2) { return (T1 + T2) / 2.0; }
};

bool MinNonZero(float i, float j) {
  float min_allowable = 0.1;
  if (j < min_allowable) {
    return 1;
  } else if (i < min_allowable) {
    return 0;
  } else {
    return i < j;
  }
  return i < j;
}

void ScaleVector(vector<float> &vec, float factor) {
  for (auto &element: vec) {
    element *= factor;
  }
}

void DumpVector(const vector<float> &vec) {
  for (auto element: vec) cout << element << endl;
}

string TruncateFloatToString(float value, unsigned n_decimal_places) {
  std::stringstream stream;
  stream << std::fixed << std::setprecision(n_decimal_places) << value;
  return stream.str();
}

}

int Plotter::PlotLumiCurrent(const vector<float> &lumi_arg,
                             const vector<float> &current_arg,
                             string run_name,
                             string channel_name,
                             const LumiCurrentPlotOptions &plot_options,
                             string output_dir,
                             FitResults &fit_results) {
  TCanvas canvas;
  canvas.cd();

  vector<float> lumi;
  vector<float> current;
  auto num_points_arg = lumi_arg.size();
  float epsilon = 0.00000000001;

  for (unsigned i=0; i < num_points_arg; ++i) {
    if (lumi_arg.at(i) > epsilon && current_arg.at(i) > epsilon &&
        !(run_name == "203934" && i == num_points_arg - 1)) {
      lumi.push_back(lumi_arg.at(i));
      current.push_back(current_arg.at(i));
    }
  }

  auto num_points = lumi.size();

  if (abs(plot_options.x_scale - 1) < epsilon) {
    ScaleVector(lumi, plot_options.x_scale);
  }
  if (abs(plot_options.y_scale - 1) < epsilon) {
    ScaleVector(current, plot_options.y_scale);
  }

  Float_t *lumi_arr = new Float_t[num_points];
  Float_t *current_arr = new Float_t[num_points];
  for (unsigned iPoint=0; iPoint < num_points; ++iPoint) {
    lumi_arr[iPoint] = lumi.at(iPoint);
    current_arr[iPoint] = current.at(iPoint);
  }

  TGraph graph(num_points, lumi_arr, current_arr);
  //graph.Draw(plot_options.draw_options.c_str());

  graph.GetXaxis()->SetTitle(plot_options.x_title.c_str());
  graph.GetYaxis()->SetTitle(plot_options.y_title.c_str());

  /*
  float x_min;
  float x_max;
  float y_min;
  float y_max;
  */

  // SetLimits doesn't work properly for Y axis, nor SetRangeUser for X axis
  //   oooookay, ROOT
  if (plot_options.x_auto_range) {
    auto lumi_max = *(std::max_element(lumi.begin(), lumi.end()));
    auto lumi_min_nonzero = *(std::min_element(lumi.begin(), lumi.end()));//, MinNonZero));
    graph.GetXaxis()->SetRangeUser(0.8*lumi_min_nonzero, 1.2*lumi_max);
  } else {
    graph.GetXaxis()->SetRangeUser(plot_options.x_min, plot_options.x_max);
  }
  if (plot_options.y_auto_range) {
    auto current_max = *(std::max_element(current.begin(), current.end()));
    auto current_min_nonzero = *(std::min_element(current.begin(), current.end()));//, MinNonZero));
    graph.GetYaxis()->SetRangeUser(0.8*current_min_nonzero, 1.2*current_max);
  } else {
    graph.GetYaxis()->SetRangeUser(plot_options.y_min, plot_options.y_max);
  }

  graph.SetMarkerColor(plot_options.marker_color);
  graph.SetMarkerStyle(plot_options.marker_style);
  graph.SetMarkerSize(plot_options.marker_size);

  string graph_title = "Run "+run_name+", Channel "+channel_name;
  graph.SetTitle( graph_title.c_str() );

  TPaveText fit_legend;
  if (plot_options.do_fit) {
    graph.Fit("pol1", plot_options.fit_options.c_str());
    TF1 *fit = graph.GetFunction("pol1");
    fit_results.slope = fit->GetParameter(1);
    fit_results.slope_err = fit->GetParError(1);
    fit_results.intercept = fit->GetParameter(0);
    fit_results.intercept_err = fit->GetParError(0);
    fit_results.chi_squared = fit->GetChisquare();
    fit_results.nDoF = fit->GetNDF();
    // Convert to calibration used to get lumi = calib*current
    //   current*y_scale = m*lumi*x_scale + intercept
    //   lumi = 1 / (x_scale*m) * (current*y_scale + intercept)
    fit_results.calibration_slope = plot_options.y_scale /
                                   (fit_results.slope*plot_options.x_scale);
    fit_results.calibration_intercept = (-1)*fit_results.intercept /
                                       (fit_results.slope*plot_options.x_scale);

    fit->SetLineColor(plot_options.fit_line_color);
    fit->SetLineWidth(plot_options.fit_line_width);

    if (plot_options.fit_show_legend) {
      float goodness_of_fit = fit_results.chi_squared / fit_results.nDoF;
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
  }

  graph.Draw(plot_options.draw_options.c_str());
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

int Plotter::SaveFitResultsToTree(const std::map<string, FitResults> &fit_results,
                                  string run_name,
                                  string output_dir) {
// Saves FCal current vs. BCM lumi fit parameters to a root file for a given
//   run. These parameters are used to derive the FCal lumi calibration.

  string write_dir = output_dir+"fit_results/";
  int err = system( ("mkdir -p "+write_dir).c_str() );

  TFile file((write_dir+run_name+".root").c_str(), "recreate");
  TTree tree;

  string channel_name;
  float slope;
  float slope_err;
  float intercept;
  float intercept_err;
  float chi_squared;
  int nDoF;
  bool is_short;

  tree.Branch("channel_name", &channel_name);
  tree.Branch("slope", &slope);
  tree.Branch("slope_err", &slope_err);
  tree.Branch("intercept", &intercept);
  tree.Branch("intercept_err", &intercept_err);
  tree.Branch("chi_squared", &chi_squared);
  tree.Branch("nDoF", &nDoF);
  tree.Branch("is_short", &is_short);

  for (const auto& this_channel_results: fit_results) {
    const auto& result = this_channel_results.second;

    channel_name = this_channel_results.first;
    slope = result.slope;
    slope_err = result.slope_err;
    intercept = result.intercept;
    intercept_err = result.intercept_err;
    chi_squared = result.chi_squared;
    nDoF = result.nDoF;
    is_short = result.is_short;

    tree.Fill();
  }

  tree.Write();

  return 0;
}

int Plotter::SaveCalibrationToText(const std::map<string, FitResults> &fit_results,
                                   string run_name,
                                   string output_dir) {
// Saves the FCal lumi calibration data space-delimited in a plaintext file:
//     
//     ChannelName slope intercept

  string write_dir = output_dir+"fit_results/calibrations/";
  string out_filepath = write_dir + "calib_" + run_name + ".dat";
  int err = system( ("mkdir -p "+write_dir).c_str() );

  std::ofstream out_file(out_filepath);
  if (!out_file.is_open()) {
    cerr << "ERROR: in SaveCalibrationToText(...):" << endl;
    cerr << "\tcould not open file \'" << out_filepath << "\'" << endl;
    return 1;
  }

  for (const auto& this_channel_results: fit_results) {
    const auto& result = this_channel_results.second;

    string channel_name = this_channel_results.first;
    float slope = result.calibration_slope;
    float intercept = result.calibration_intercept;

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

  TLegend legend(0.49,0.65,0.85,0.85);
  legend.SetFillColor(0);
  legend.SetBorderSize(0);
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
  if (plot_options.x_auto_range) {
    low_bin = runs_data.begin()->second.timestamp();
    high_bin = runs_data.end()->second.timestamp();
    n_bins = 4320;
  } else {
    low_bin = plot_options.low_bin;
    high_bin = plot_options.high_bin;
    n_bins = (high_bin - low_bin) / plot_options.bin_size;
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
  stack.Draw(plot_options.draw_options.c_str());
  FormatTHStack(plot_options, stack);

  legend.Draw();

  canvas.Update();
  canvas.Print( (write_dir+"mu_stability.pdf").c_str() );

  return 0;
}
