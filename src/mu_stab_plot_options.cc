#include <string>

#include "Rtypes.h"

#include "json_reader.h"
#include "mu_stab_plot_options.h"

using std::string;

MuStabPlotOptions::MuStabPlotOptions(string params_filepath) {
  JSONReader parameter_file(params_filepath);
  string node_name = "plot_options.mu_stability.";

  base_output_dir_ = parameter_file.get<string>("output_dirs.base") +
                     parameter_file.get<string>("output_dirs.mu_stability");

  rootfiles_output_dir_ = base_output_dir_ +
                          parameter_file.get<string>(
                              node_name+"output_dirs.rootfiles");

  draw_options_ = parameter_file.get<string>(node_name+"draw_options");
  detector_name_ = parameter_file.get<string>("detector");

  marker_color_A_ = parameter_file.get<int>(node_name+"marker.A.color");
  marker_size_A_ = parameter_file.get<Float_t>(node_name+"marker.A.size");
  marker_style_A_ = parameter_file.get<int>(node_name+"marker.A.style");
  marker_color_C_ = parameter_file.get<int>(node_name+"marker.C.color");
  marker_size_C_ = parameter_file.get<Float_t>(node_name+"marker.C.size");
  marker_style_C_ = parameter_file.get<int>(node_name+"marker.C.style");
  marker_color_Avg_ = parameter_file.get<int>(node_name+"marker.Avg.color");
  marker_size_Avg_ = parameter_file.get<Float_t>(node_name+"marker.Avg.size");
  marker_style_Avg_ = parameter_file.get<int>(node_name+"marker.Avg.style");

  x_auto_range_ = parameter_file.get<bool>(node_name+"x.auto_range");
  low_bin_ = parameter_file.get<UInt_t>(node_name+"x.min");
  high_bin_ = parameter_file.get<UInt_t>(node_name+"x.max");
  n_bins_ = parameter_file.get<UInt_t>(node_name+"x.n_bins");
  x_title_ = parameter_file.get<string>(node_name+"x.title");

  y_auto_range_ = parameter_file.get<bool>(node_name+"y.auto_range");
  y_min_ = parameter_file.get<Float_t>(node_name+"y.min");
  y_max_ = parameter_file.get<Float_t>(node_name+"y.max");
  y_title_ = parameter_file.get<string>(node_name+"y.title");
}
