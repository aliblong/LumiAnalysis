#include "scatter_plot_options.h"

#include <string>

#include "Rtypes.h"

#include "json_reader.h"

using std::string;

ScatterPlotOptions::ScatterPlotOptions(const JSONReader& params, const std::string& node, bool verbose) :
  x_(params, node + "x.", verbose),
  y_(params, node + "y.", verbose)
{
  output_dir_ = params.get<string>("base_output_dir", "output/", verbose) + params.get<string>(node+"output_dirs.base", "scatter_default/", verbose) + params.get<string>(node+"output_dirs.plots", "plots/", verbose);

  draw_options_ = params.get<string>(node+"draw_options", "", verbose);
  title_ = params.get<string>(node+"title", "", verbose);

  marker_colors_ = params.get_vector<Int_t>(node+"marker.colors", {2}, verbose);
  marker_sizes_ = params.get_vector<Float_t>(node+"marker.sizes", {0.3}, verbose);
  marker_styles_ = params.get_vector<Int_t>(node+"marker.styles", {20}, verbose);

}

AxisOptions::AxisOptions(const JSONReader& params, const std::string& node, bool verbose)
{
  scale_ = params.get<Float_t>(node+"scale", 1.0, verbose);
  rel_error_ = params.get<Float_t>(node+"rel_error", 0.1, verbose);
  auto_range_ = params.get<bool>(node+"auto_range", true, verbose);
  min_ = params.get<Float_t>(node+"min", 0.0, verbose);
  max_ = params.get<Float_t>(node+"max", 1.0, verbose);
  title_ = params.get<string>(node+"title", "", verbose);
}
