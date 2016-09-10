#include <string>

#include "Rtypes.h"

#include "json_reader.h"
#include "lumi_current_plot_options.h"

using std::string;

// Decline to handle errors from missing node in both input and default parameter trees
LumiCurrentPlotOptions::LumiCurrentPlotOptions(const JSONReader& params, const string& node) :
  ScatterPlotOptions(params, node)
{
  auto base_output_dir = *params.get<string>("base_output_dir") +
                         *params.get<string>(node+"output_dirs.base");

  plots_dir_ = base_output_dir +
               *params.get<string>(node+"output_dirs.plots");

  auto LC_fit_output_node = node + "output_dirs.fit_results.";
  auto base_fit_results_dir = base_output_dir +
                              *params.get<string>(
                                LC_fit_output_node+"base");
  raw_fit_results_dir_ = base_fit_results_dir +
                         *params.get<string>(
                           LC_fit_output_node+"raw");
  calibration_results_dir_ = base_fit_results_dir +
                             *params.get<string>(
                               LC_fit_output_node+"calibrations");
  geometric_results_dir_ = base_fit_results_dir +
                           *params.get<string>(
                             LC_fit_output_node+"geometric");

  do_individual_ = *params.get<bool>(node+"do_individual");
  do_sum_ = *params.get<bool>(node+"do_sum");
  do_fit_ = *params.get<bool>(node+"fit.do");
  print_plots_ = *params.get<bool>(node+"print_plots");

  // Fit options
  if (do_fit_) {
    string fit_verbosity_option;
    int fit_verbosity = *params.get<int>(node+"fit.verbose");
    if (fit_verbosity == 2) {
      fit_verbosity_option = "V";
    }
    else if (fit_verbosity == 1) {
      fit_verbosity_option = "";
    }
    else { // verbosity 0
      fit_verbosity_option = "Q";
    }
    string fit_show_option;
    if (*params.get<bool>(node+"fit.show")) {
      fit_show_option = "";
    }
    else {
      fit_show_option = "0";
    }
    string additional_options = *params.get<string>(node+"fit.options");

    fit_options_ = fit_verbosity_option + fit_show_option + additional_options;

    fit_fix_intercept_ = *params.get<bool>(node+"fit.fix_intercept");

    fit_show_legend_ = *params.get<bool>(node+"fit.show_legend");
    fit_line_color_ = *params.get<int>(node+"fit.line_color");
    fit_line_width_ = *params.get<Float_t>(node+"fit.line_width");
  }
}
