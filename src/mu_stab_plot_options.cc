#include <string>

#include "Rtypes.h"

#include "json_reader.h"
#include "mu_stab_plot_options.h"

using std::string;

MuStabPlotOptions::MuStabPlotOptions(string params_filepath) {
  JSONReader parameter_file(params_filepath);
  string node_name = "plot_options.mu_stability.";

  draw_options = parameter_file.get<string>(node_name+"draw_options");

  marker_color_A = parameter_file.get<int>(node_name+"marker.A.color");
  marker_size_A = parameter_file.get<float>(node_name+"marker.A.size");
  marker_style_A = parameter_file.get<int>(node_name+"marker.A.style");
  marker_color_C = parameter_file.get<int>(node_name+"marker.C.color");
  marker_size_C = parameter_file.get<float>(node_name+"marker.C.size");
  marker_style_C = parameter_file.get<int>(node_name+"marker.C.style");
  marker_color_Avg = parameter_file.get<int>(node_name+"marker.Avg.color");
  marker_size_Avg = parameter_file.get<float>(node_name+"marker.Avg.size");
  marker_style_Avg = parameter_file.get<int>(node_name+"marker.Avg.style");

  x_auto_range = parameter_file.get<bool>(node_name+"x.auto_range");
  low_bin = parameter_file.get<UInt_t>(node_name+"x.min");
  high_bin = parameter_file.get<UInt_t>(node_name+"x.max");
  bin_size = parameter_file.get<UInt_t>(node_name+"x.bin_size");
  x_title = parameter_file.get<string>(node_name+"x.title");

  y_auto_range = parameter_file.get<bool>(node_name+"y.auto_range");
  y_min = parameter_file.get<float>(node_name+"y.min");
  y_max = parameter_file.get<float>(node_name+"y.max");
  y_title = parameter_file.get<string>(node_name+"y.title");
}
