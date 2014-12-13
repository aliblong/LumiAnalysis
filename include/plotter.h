#ifndef LUMIANALYSIS_INCLUDE_PLOTTER_H_
#define LUMIANALYSIS_INCLUDE_PLOTTER_H_

#include <string>
#include <vector>
#include <map>

#include "fit_results.h"
#include "lumi_current_plot_options.h"
#include "mu_stab_plot_options.h"
#include "single_run_data.h"

namespace Plotter {
  int GeometricAnalysisOfFitResults(
                       const std::map<std::string, FitResults> &fit_results,
                       std::string run_name,
                       std::string output_dir);

  int PlotLumiCurrent(const std::vector<Float_t> &lumi,
                      const std::vector<Float_t> &current,
                      std::string run_name,
                      std::string channel_name,
                      const LumiCurrentPlotOptions &plot_options,
                      std::string output_dir,
                      FitResults &fit_results);

  int PlotMuStability(const std::map<std::string, SingleRunData> &runs_data,
                      const MuStabPlotOptions &plot_options,
                      std::string output_dir);

  int WriteFitResultsToTree(const std::map<std::string, FitResults> &fit_results,
                            std::string run_name,
                            std::string output_name);

  int WriteCalibrationToText(const std::map<std::string, FitResults> &fit_results,
                             std::string run_name,
                             std::string output_dir);
}

#endif
