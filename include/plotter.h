#ifndef LUMIANALYSIS_INCLUDE_PLOTTER_H_
#define LUMIANALYSIS_INCLUDE_PLOTTER_H_

#include <string>
#include <vector>

using std::string;
using std::vector;

class Plotter {
 public:
  Plotter(){};
  ~Plotter(){};
  int PlotLumiCurrent(vector<float> lumi,
                      vector<float> current,
                      string run_name,
                      string channel_name,
                      string output_dir);
};

#endif
