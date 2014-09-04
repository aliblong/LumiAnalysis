#ifndef LUMIANALYSIS_INCLUDE_PLOTTER_H_
#define LUMIANALYSIS_INCLUDE_PLOTTER_H_

#include <string>
#include <vector>
#include <map>

#include "plot_options.h"
#include "fit_results.h"

class Plotter {
  typedef typename std::string string;
 public:
  Plotter(){};
  ~Plotter(){};

  int PlotLumiCurrent(const std::vector<float> &lumi,
                      const std::vector<float> &current,
                      string run_name,
                      string channel_name,
                      string output_dir,
                      const PlotOptions &plot_options,
                      FitResults &fit_results);

  int SaveFitResults(const std::map<string, FitResults> &fit_results,
                     string run_name,
                     string output_name);
};

#endif
