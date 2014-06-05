#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <memory>
#include <cstdlib>
#include <algorithm>
#include <iomanip>

#include "TAxis.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TGraph.h"
#include "TPaveText.h"

#include "plotter.h"
#include "plot_options.h"

using std::cout;
using std::cerr;
using std::endl;

using std::string;
using std::vector;

namespace {

bool MinNonZero(float i, float j) {
  float min_allowable = 0.1;
  if (j < min_allowable) {
    return 1;
  } else if (i < min_allowable) {
    return 0;
  } else {
    return i < j;
  }
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
                             string output_dir,
                             const PlotOptions &plot_options) {
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
    auto lumi_min_nonzero = *(std::min_element(lumi.begin(), lumi.end(), MinNonZero));
    graph.GetXaxis()->SetRangeUser(0.8*lumi_min_nonzero, 1.2*lumi_max);
  } else {
    graph.GetXaxis()->SetLimits(plot_options.x_min, plot_options.x_max);
  }
  if (plot_options.y_auto_range) {
    auto current_max = *(std::max_element(current.begin(), current.end()));
    auto current_min_nonzero = *(std::min_element(current.begin(), current.end(), MinNonZero));
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
    fit->SetLineColor(plot_options.fit_line_color);
    fit->SetLineWidth(plot_options.fit_line_width);

    if (plot_options.fit_show_legend) {
      string chi2 = TruncateFloatToString(fit->GetChisquare(), 2);
      string slope = TruncateFloatToString(fit->GetParameter(1), 2);
      string slope_err = TruncateFloatToString(fit->GetParError(1), 2);
      string intercept = TruncateFloatToString(fit->GetParameter(0), 2);
      string intercept_err = TruncateFloatToString(fit->GetParError(0), 2);
      fit_legend.AddText( ("#chi^{2} = "+chi2).c_str() );
      fit_legend.AddText( ("slope = "+slope+" #pm "+slope_err).c_str() );
      fit_legend.AddText( ("constant = "+intercept+" #pm "+intercept_err).c_str() );
      fit_legend.SetX1NDC(0.15);
      fit_legend.SetX2NDC(0.4);
      fit_legend.SetY1NDC(0.65);
      fit_legend.SetY2NDC(0.9);
      fit_legend.SetTextSize(0.035);
      fit_legend.SetOption("");
    }
  }

  graph.Draw(plot_options.draw_options.c_str());
  fit_legend.Draw();

  canvas.Update();
  string write_dir = output_dir+"lumi_current/"+run_name+"/";
  int err_system = system( ("mkdir -p "+write_dir).c_str() );

  canvas.Print( (write_dir+channel_name+".png").c_str() );
  TFile *this_file = TFile::Open("graph_file.root", "RECREATE");
  graph.Write();
  this_file->Close();

  delete [] lumi_arr;
  delete [] current_arr;
  return 0;
}
