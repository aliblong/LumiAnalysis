#include <string>

#include "Rtypes.h"

#include "json_reader.h"
#include "lumi_current_plot_options.h"

using std::string;

LumiCurrentPlotOptions::LumiCurrentPlotOptions(const JSONReader& params) {
  string LC_options_node = "plot_options.lumi_current.";

  auto base_output_dir = params.get<string>("output_dirs.base") +
                         params.get<string>("output_dirs.lumi_current");

  plots_dir_ = base_output_dir +
               params.get<string>(LC_options_node+"output_dirs.plots");

  auto LC_fit_output_node = LC_options_node + "output_dirs.fit_results.";
  auto base_fit_results_dir = base_output_dir +
                              params.get<string>(
                                LC_fit_output_node+"base");
  raw_fit_results_dir_ = base_fit_results_dir +
                         params.get<string>(
                           LC_fit_output_node+"raw");
  calibration_results_dir_ = base_fit_results_dir +
                             params.get<string>(
                               LC_fit_output_node+"calibrations");
  geometric_results_dir_ = base_fit_results_dir +
                           params.get<string>(
                             LC_fit_output_node+"geometric");

  do_individual_ = params.get<bool>(LC_options_node+"do_individual");
  do_sum_ = params.get<bool>(LC_options_node+"do_sum");
  do_fit_ = params.get<bool>(LC_options_node+"fit.do");
  print_plots_ = params.get<bool>(LC_options_node+"print_plots");

  // Fit options
  string fit_verbosity_option;
  int fit_verbosity = params.get<int>(LC_options_node+"fit.verbose");
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
  if (params.get<bool>(LC_options_node+"fit.show")) {
    fit_show_option = "";
  }
  else {
    fit_show_option = "0";
  }
  string additional_options = params.get<string>(LC_options_node+"fit.options");

  fit_options_ = fit_verbosity_option + fit_show_option + additional_options;

  fit_fix_intercept_ = params.get<bool>(LC_options_node+"fit.fix_intercept");

  fit_show_legend_ = params.get<bool>(LC_options_node+"fit.show_legend");
  fit_line_color_ = params.get<int>(LC_options_node+"fit.line_color");
  fit_line_width_ = params.get<Float_t>(LC_options_node+"fit.line_width");

  draw_options_ = params.get<string>(LC_options_node+"draw_options");

  marker_color_ = params.get<int>(LC_options_node+"marker.color");
  marker_style_ = params.get<int>(LC_options_node+"marker.style");
  marker_size_ = params.get<Float_t>(LC_options_node+"marker.size");

  x_scale_ = params.get<Float_t>(LC_options_node+"x.scale");
  x_rel_error_ = params.get<Float_t>(LC_options_node+"x.rel_error");
  x_auto_range_ = params.get<bool>(LC_options_node+"x.auto_range");
  x_min_ = params.get<Float_t>(LC_options_node+"x.min");
  x_max_ = params.get<Float_t>(LC_options_node+"x.max");
  x_title_ = params.get<string>(LC_options_node+"x.title");

  y_scale_ = params.get<Float_t>(LC_options_node+"y.scale");
  y_rel_error_ = params.get<Float_t>(LC_options_node+"y.rel_error");
  y_auto_range_ = params.get<bool>(LC_options_node+"y.auto_range");
  y_min_ = params.get<Float_t>(LC_options_node+"y.min");
  y_max_ = params.get<Float_t>(LC_options_node+"y.max");
  y_title_ = params.get<string>(LC_options_node+"y.title");
}
