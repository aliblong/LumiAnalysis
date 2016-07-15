#include "scatter_plot_options.h"

#include <string>

#include "Rtypes.h"

#include "json_reader.h"

using std::string;

ScatterPlotOptions::ScatterPlotOptions(const JSONReader& params, const std::string& node) :
  x_(params, node + "x."),
  y_(params, node + "y.")
{
  output_dir_ = params.get<string>("base_output_dir") + params.get<string>(node+"output_dirs.base") + params.get<string>(node+"output_dirs.plots");

  draw_options_ = params.get<string>(node+"draw_options");
  title_ = params.get<string>(node+"title", "");

  marker_colors_ = params.get_vector<Int_t>(node+"marker.colors");
  marker_sizes_ = params.get_vector<Float_t>(node+"marker.sizes");
  marker_styles_ = params.get_vector<Int_t>(node+"marker.styles");

}

AxisOptions::AxisOptions(const JSONReader& params, const std::string& node)
{
  scale_ = params.get<Float_t>(node+"scale");
  rel_error_ = params.get<Float_t>(node+"rel_error");
  auto_range_ = params.get<bool>(node+"auto_range");
  min_ = params.get<Float_t>(node+"min");
  max_ = params.get<Float_t>(node+"max");
  title_ = params.get<string>(node+"title");
}
