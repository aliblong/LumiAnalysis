#include <string>

#include "Rtypes.h"

#include "json_reader.h"
#include "mu_stab_plot_options.h"

using std::string;

MuStabPlotOptions::MuStabPlotOptions(const JSONReader& params, const string& node, bool verbose) {
  base_output_dir_ = params.get<string>("output_dirs.base", "output/", verbose) +
                     params.get<string>("output_dirs.mu_stability", "mu_stab/", verbose);

  rootfiles_output_dir_ = base_output_dir_ +
                          params.get<string>(
                              node+"output_dirs.rootfiles", "rootfiles/", verbose);

  draw_options_ = params.get<string>(node+"draw_options", "", verbose);
  detector_name_ = params.get<string>("detector", "FCal", verbose);

  marker_color_A_ = params.get<Int_t>(node+"marker.A.color", 797, verbose);
  marker_size_A_ = params.get<Float_t>(node+"marker.A.size", 0.9, verbose);
  marker_style_A_ = params.get<Int_t>(node+"marker.A.style", 22, verbose);
  marker_color_C_ = params.get<Int_t>(node+"marker.C.color", 9, verbose);
  marker_size_C_ = params.get<Float_t>(node+"marker.C.size", 0.9, verbose);
  marker_style_C_ = params.get<Int_t>(node+"marker.C.style", 22, verbose);
  marker_color_Avg_ = params.get<Int_t>(node+"marker.Avg.color", 20, verbose);
  marker_size_Avg_ = params.get<Float_t>(node+"marker.Avg.size", 0.9, verbose);
  marker_style_Avg_ = params.get<Int_t>(node+"marker.Avg.style", 22, verbose);

  x_auto_range_ = params.get<bool>(node+"x.auto_range", true, verbose);
  low_bin_ = params.get<UInt_t>(node+"x.min", 0, verbose);
  high_bin_ = params.get<UInt_t>(node+"x.max", 1, verbose);
  n_bins_ = params.get<UInt_t>(node+"x.n_bins", 1, verbose);
  x_title_ = params.get<string>(node+"x.title", "default", verbose);

  y_auto_range_ = params.get<bool>(node+"y.auto_range", true, verbose);
  y_min_ = params.get<Float_t>(node+"y.min", 0.0, verbose);
  y_max_ = params.get<Float_t>(node+"y.max", 1.0, verbose);
  y_title_ = params.get<string>(node+"y.title", "default", verbose);
}
