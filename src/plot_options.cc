#include <string>

#include "json_reader.h"
#include "plot_options.h"

using std::string;

PlotOptions::PlotOptions(string params_filepath, string plot_type) {
  JSONReader parameter_file(params_filepath);
  string node_name = "plot_options."+plot_type+".";

  do_fit = parameter_file.get<bool>(node_name+"fit.do");

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
  fit_options = fit_verbosity_option + fit_show_option;

  fit_show_legend = parameter_file.get<bool>(node_name+"fit.show_legend");
  fit_line_color = parameter_file.get<int>(node_name+"fit.line_color");
  fit_line_width = parameter_file.get<float>(node_name+"fit.line_width");

  draw_options = parameter_file.get<string>(node_name+"draw_options");

  marker_color = parameter_file.get<int>(node_name+"marker.color");
  marker_style = parameter_file.get<int>(node_name+"marker.style");
  marker_size = parameter_file.get<float>(node_name+"marker.size");

  x_scale = parameter_file.get<float>(node_name+"x.scale");
  x_auto_range = parameter_file.get<bool>(node_name+"x.auto_range");
  x_min = parameter_file.get<float>(node_name+"x.min");
  x_max = parameter_file.get<float>(node_name+"x.max");
  x_title = parameter_file.get<string>(node_name+"x.title");

  y_scale = parameter_file.get<float>(node_name+"y.scale");
  y_auto_range = parameter_file.get<bool>(node_name+"y.auto_range");
  y_min = parameter_file.get<float>(node_name+"y.min");
  y_max = parameter_file.get<float>(node_name+"y.max");
  y_title = parameter_file.get<string>(node_name+"y.title");
}
