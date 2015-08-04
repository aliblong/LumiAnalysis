#include <string>

#include "Rtypes.h"

#include "json_reader.h"
#include "mu_dep_plot_options.h"

using std::string;

MuDepPlotOptions::MuDepPlotOptions(const string& params_filepath) {
  JSONReader parameter_file(params_filepath);

  string LC_options_node = "plot_options.mu_dependence.";

  auto base_output_dir = parameter_file.get<string>("output_dirs.base") +
                         parameter_file.get<string>("output_dirs.mu_dependence");

  plots_dir_ = base_output_dir +
               parameter_file.get<string>(LC_options_node+"output_dirs.plots");

  draw_options_ = parameter_file.get<string>(LC_options_node+"draw_options");

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
