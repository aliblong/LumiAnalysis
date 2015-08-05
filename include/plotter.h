#ifndef LUMIANALYSIS_INCLUDE_PLOTTER_H_
#define LUMIANALYSIS_INCLUDE_PLOTTER_H_

#include <string>
#include <vector>

#include "TGraphErrors.h"

#include "boost/container/flat_map.hpp"
#include "boost/expected/expected.hpp"

#include "error.h"
#include "fit_results.h"
#include "lumi_current_plot_options.h"
#include "mu_lumi_dep_plot_options.h"
#include "mu_stab_plot_options.h"
#include "single_run_data.h"
#include "void.h"

namespace Plotter {
  Error::Expected<Void> GeometricAnalysisOfFitResults(
      const boost::container::flat_map<std::string, FitResults>& fit_results,
      const LumiCurrentPlotOptions& options);

  Error::Expected<FitResults> PlotLumiCurrent(
      const std::vector<std::array<Float_t, 2>>& points,
      const LumiCurrentPlotOptions& options);

  Error::Expected<Void> PlotMuLumiDependence(
      const std::vector<std::array<Float_t, 2>>& points,
      const MuLumiDepPlotOptions& options);

  Error::Expected<Void> PlotMuStability(
      const boost::container::flat_map<std::string, SingleRunData>& runs_data,
      const MuStabPlotOptions& plot_options);

  Error::Expected<Void> WriteFitResultsToTree(
      const boost::container::flat_map<std::string, FitResults>& fit_results,
      const LumiCurrentPlotOptions& options);

  Error::Expected<Void> WriteCalibrationToText(
      const boost::container::flat_map<std::string, FitResults>& fit_results,
      const LumiCurrentPlotOptions& options);
}

/*
  Error::Expected<Void> PlotLumiTotalCurrent(const std::vector<Float_t>& lumi,
                           const std::vector<Float_t>& current_A,
                           const std::vector<Float_t>& current_C,
                           std::string run_name,
                           const LumiCurrentPlotOptions& plot_options,
                           std::string output_dir);
                           */

#endif
