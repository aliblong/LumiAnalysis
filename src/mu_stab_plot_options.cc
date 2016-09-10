#include <string>

#include "Rtypes.h"

#include "json_reader.h"
#include "mu_stab_plot_options.h"

using std::string;

MuStabPlotOptions::MuStabPlotOptions(const JSONReader& params, const string& node) {
  base_output_dir_ = *params.get<string>("output_dirs.base") +
                     *params.get<string>("output_dirs.mu_stability");

  rootfiles_output_dir_ = base_output_dir_ +
                          *params.get<string>(
                              node+"output_dirs.rootfiles");

  draw_options_ = *params.get<string>(node+"draw_options");
  detector_name_ = *params.get<string>("detector");

  marker_color_A_ = *params.get<Int_t>(node+"marker.A.color");
  marker_size_A_ = *params.get<Float_t>(node+"marker.A.size");
  marker_style_A_ = *params.get<Int_t>(node+"marker.A.style");
  marker_color_C_ = *params.get<Int_t>(node+"marker.C.color");
  marker_size_C_ = *params.get<Float_t>(node+"marker.C.size");
  marker_style_C_ = *params.get<Int_t>(node+"marker.C.style");
  marker_color_Avg_ = *params.get<Int_t>(node+"marker.Avg.color");
  marker_size_Avg_ = *params.get<Float_t>(node+"marker.Avg.size");
  marker_style_Avg_ = *params.get<Int_t>(node+"marker.Avg.style");

  x_auto_range_ = *params.get<bool>(node+"x.auto_range");
  low_bin_ = *params.get<UInt_t>(node+"x.min");
  high_bin_ = *params.get<UInt_t>(node+"x.max");
  n_bins_ = *params.get<UInt_t>(node+"x.n_bins");
  x_title_ = *params.get<string>(node+"x.title");

  y_auto_range_ = *params.get<bool>(node+"y.auto_range");
  y_min_ = *params.get<Float_t>(node+"y.min");
  y_max_ = *params.get<Float_t>(node+"y.max");
  y_title_ = *params.get<string>(node+"y.title");
}
