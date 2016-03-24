#include "mu_lumi_dep_plot_options.h"

#include <string>

#include "Rtypes.h"

#include "json_reader.h"

using std::string;

MuLumiDepPlotOptions::MuLumiDepPlotOptions(const string& params_filepath, const std::string& node) {
  JSONReader parameter_file(params_filepath);

  string LC_options_node = node;

  auto base_output_dir = parameter_file.get<string>("output_dirs.base") +
                         parameter_file.get<string>("output_dirs.mu_dependence");

  plots_dir_ = base_output_dir +
               parameter_file.get<string>(LC_options_node+"output_dirs.plots");

  // Fit options
  do_fit_ = parameter_file.get<bool>(LC_options_node+"fit.do");
  string fit_verbosity_option;
  int fit_verbosity = parameter_file.get<int>(LC_options_node+"fit.verbose");
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
  if (parameter_file.get<bool>(LC_options_node+"fit.show")) {
    fit_show_option = "";
  }
  else {
    fit_show_option = "0";
  }
  string additional_options = parameter_file.get<string>(LC_options_node+"fit.options");
  fit_options_ = fit_verbosity_option + fit_show_option + additional_options;

  fit_fix_intercept_ = parameter_file.get<bool>(LC_options_node+"fit.fix_intercept");
  fit_show_legend_ = parameter_file.get<bool>(LC_options_node+"fit.show_legend");
  fit_line_color_ = parameter_file.get<int>(LC_options_node+"fit.line_color");
  fit_line_width_ = parameter_file.get<Float_t>(LC_options_node+"fit.line_width");

  file_name_ = parameter_file.get<string>(LC_options_node+"file_name");
  draw_options_ = parameter_file.get<string>(LC_options_node+"draw_options");
  title_ = parameter_file.get<string>(LC_options_node+"title");

  marker_color_ = parameter_file.get<int>(LC_options_node+"marker.color");
  marker_style_ = parameter_file.get<int>(LC_options_node+"marker.style");
  marker_size_ = parameter_file.get<Float_t>(LC_options_node+"marker.size");

  x_scale_ = parameter_file.get<Float_t>(LC_options_node+"x.scale");
  x_rel_error_ = parameter_file.get<Float_t>(LC_options_node+"x.rel_error");
  x_auto_range_ = parameter_file.get<bool>(LC_options_node+"x.auto_range");
  x_min_ = parameter_file.get<Float_t>(LC_options_node+"x.min");
  x_max_ = parameter_file.get<Float_t>(LC_options_node+"x.max");
  x_title_ = parameter_file.get<string>(LC_options_node+"x.title");

  y_scale_ = parameter_file.get<Float_t>(LC_options_node+"y.scale");
  y_rel_error_ = parameter_file.get<Float_t>(LC_options_node+"y.rel_error");
  y_auto_range_ = parameter_file.get<bool>(LC_options_node+"y.auto_range");
  y_min_ = parameter_file.get<Float_t>(LC_options_node+"y.min");
  y_max_ = parameter_file.get<Float_t>(LC_options_node+"y.max");
  y_title_ = parameter_file.get<string>(LC_options_node+"y.title");
}
