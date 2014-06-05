#ifndef LUMIANALYSIS_INCLUDE_PLOTTER_H_
#define LUMIANALYSIS_INCLUDE_PLOTTER_H_

#include <string>
#include <vector>

#include "plot_options.h"

using std::string;
using std::vector;

class Plotter {
 public:
  Plotter(){};
  ~Plotter(){};
  int PlotLumiCurrent(const vector<float> &lumi,
                      const vector<float> &current,
                      string run_name,
                      string channel_name,
                      string output_dir,
                      const PlotOptions &plot_options);
};

#endif
