#include <string>

#include "Rtypes.h"

#include "json_reader.h"
#include "lumi_current_plot_options.h"

using std::string;

LumiCurrentPlotOptions::LumiCurrentPlotOptions(string params_filepath) {
  JSONReader parameter_file(params_filepath);
  string node_name = "plot_options.lumi_current.";

  do_individual_ = parameter_file.get<bool>(node_name+"do_individual");
  do_sum_ = parameter_file.get<bool>(node_name+"do_sum");
  do_fit_ = parameter_file.get<bool>(node_name+"fit.do");

  // Fit options
  string fit_verbosity_option;
  int fit_verbosity = parameter_file.get<int>(node_name+"fit.verbose");
  if (fit_verbosity == 2) {
    fit_verbosity_option = "V";
  } else if (fit_verbosity == 1) {
    fit_verbosity_option = "";
  } else { // verbosity 0
    fit_verbosity_option = "Q";
  }
  string fit_show_option;
  if (parameter_file.get<bool>(node_name+"fit.show")) {
    fit_show_option = "";
  } else {
    fit_show_option = "0";
  }
  string additional_options = parameter_file.get<string>(node_name+"fit.options");

  fit_options_ = fit_verbosity_option + fit_show_option + additional_options;

  fit_fix_intercept_ = parameter_file.get<bool>(node_name+"fit.fix_intercept");

  fit_show_legend_ = parameter_file.get<bool>(node_name+"fit.show_legend");
  fit_line_color_ = parameter_file.get<int>(node_name+"fit.line_color");
  fit_line_width_ = parameter_file.get<Float_t>(node_name+"fit.line_width");

  draw_options_ = parameter_file.get<string>(node_name+"draw_options");

  marker_color_ = parameter_file.get<int>(node_name+"marker.color");
  marker_style_ = parameter_file.get<int>(node_name+"marker.style");
  marker_size_ = parameter_file.get<Float_t>(node_name+"marker.size");

  x_scale_ = parameter_file.get<Float_t>(node_name+"x.scale");
  x_rel_error_ = parameter_file.get<Float_t>(node_name+"x.rel_error");
  x_auto_range_ = parameter_file.get<bool>(node_name+"x.auto_range");
  x_min_ = parameter_file.get<Float_t>(node_name+"x.min");
  x_max_ = parameter_file.get<Float_t>(node_name+"x.max");
  x_title_ = parameter_file.get<string>(node_name+"x.title");

  y_scale_ = parameter_file.get<Float_t>(node_name+"y.scale");
  y_rel_error_ = parameter_file.get<Float_t>(node_name+"y.rel_error");
  y_auto_range_ = parameter_file.get<bool>(node_name+"y.auto_range");
  y_min_ = parameter_file.get<Float_t>(node_name+"y.min");
  y_max_ = parameter_file.get<Float_t>(node_name+"y.max");
  y_title_ = parameter_file.get<string>(node_name+"y.title");
}
