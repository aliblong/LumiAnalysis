#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <cstdlib>
#include <algorithm>

#include "TAxis.h"
#include "TCanvas.h"
#include "TGraph.h"

#include "plotter.h"

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

}

int Plotter::PlotLumiCurrent(vector<float> lumi,
                             vector<float> current,
                             string run_name,
                             string channel_name,
                             string output_dir) {
  TCanvas canvas;
  canvas.cd();

  auto num_points = lumi.size();

  auto current_max = *(std::max_element(current.begin(), current.end()));
  auto current_min_nonzero = *(std::min_element(current.begin(), current.end(), MinNonZero));
  auto lumi_max = *(std::max_element(lumi.begin(), lumi.end()));
  auto lumi_min_nonzero = *(std::min_element(lumi.begin(), lumi.end(), MinNonZero));

  Float_t *lumi_arr = new Float_t[num_points];
  Float_t *current_arr = new Float_t[num_points];
  for (unsigned iPoint=0; iPoint < num_points; ++iPoint) {
    lumi_arr[iPoint] = lumi.at(iPoint);
    current_arr[iPoint] = current.at(iPoint);
  }

  TGraph graph(num_points, lumi_arr, current_arr);

  graph.GetXaxis()->SetTitle("Instantaneous Luminosity [10^{30} cm^{-2} s^{-1}]");
  graph.GetXaxis()->SetRangeUser(0.8*lumi_min_nonzero, 1.2*lumi_max);
  graph.GetYaxis()->SetTitle("Current [#mu A]");
  graph.GetYaxis()->SetRangeUser(0.8*current_min_nonzero, 1.2*current_max);

  string graph_title = "Run "+run_name+", Channel "+channel_name;
  graph.SetTitle( graph_title.c_str() );
  graph.Draw("AP");

  canvas.Update();
  string write_dir = output_dir+"lumi_current/"+run_name+"/";
  int err_system = system( ("mkdir -p "+write_dir).c_str() );

  canvas.Print( (write_dir+channel_name+".png").c_str() );

  delete [] lumi_arr;
  delete [] current_arr;
  return 0;
}
