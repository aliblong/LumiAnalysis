#include <string>

#include "Rtypes.h"

#include "json_reader.h"
#include "lumi_current_plot_options.h"

using std::string;

LumiCurrentPlotOptions::LumiCurrentPlotOptions(const JSONReader& params, const string& node, bool verbose) :
  ScatterPlotOptions(params, node, verbose)
{
  auto base_output_dir = params.get<string>("base_output_dir", "output/", verbose) +
                         params.get<string>(node+"output_dirs.base", "default/", verbose);

  plots_dir_ = base_output_dir +
               params.get<string>(node+"output_dirs.plots", "plots/", verbose);

  auto LC_fit_output_node = node + "output_dirs.fit_results.";
  auto base_fit_results_dir = base_output_dir +
                              params.get<string>(
                                LC_fit_output_node+"base", "fit_results/", verbose);
  raw_fit_results_dir_ = base_fit_results_dir +
                         params.get<string>(
                           LC_fit_output_node+"raw", "raw/", verbose);
  calibration_results_dir_ = base_fit_results_dir +
                             params.get<string>(
                               LC_fit_output_node+"calibrations", "calibrations/", verbose);
  geometric_results_dir_ = base_fit_results_dir +
                           params.get<string>(
                             LC_fit_output_node+"geometric", "geometric/", verbose);

  do_individual_ = params.get<bool>(node+"do_individual", false, verbose);
  do_sum_ = params.get<bool>(node+"do_sum", false, verbose);
  do_fit_ = params.get<bool>(node+"fit.do", false, verbose);
  print_plots_ = params.get<bool>(node+"print_plots", false, verbose);

  // Fit options
  if (do_fit_) {
    string fit_verbosity_option;
    int fit_verbosity = params.get<int>(node+"fit.verbose", 2, verbose);
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
    if (params.get<bool>(node+"fit.show", false, verbose)) {
      fit_show_option = "";
    }
    else {
      fit_show_option = "0";
    }
    string additional_options = params.get<string>(node+"fit.options", "", verbose);

    fit_options_ = fit_verbosity_option + fit_show_option + additional_options;

    fit_fix_intercept_ = params.get<bool>(node+"fit.fix_intercept", false, verbose);

    fit_show_legend_ = params.get<bool>(node+"fit.show_legend", false, verbose);
    fit_line_color_ = params.get<int>(node+"fit.line_color", 1, verbose);
    fit_line_width_ = params.get<Float_t>(node+"fit.line_width", 1.0, verbose);
  }
}
