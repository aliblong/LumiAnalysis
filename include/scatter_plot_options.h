#ifndef LUMIANALYSIS_INCLUDE_SCATTER_PLOT_OPTIONS_H_
#define LUMIANALYSIS_INCLUDE_SCATTER_PLOT_OPTIONS_H_

#include <string>

#include "Rtypes.h"

#include "json_reader.h"

class AxisOptions {
  friend class ScatterPlotOptions;
 public:
  AxisOptions(const JSONReader& params, const std::string& node);
  auto scale() const { return scale_; }
  auto rel_error() const { return rel_error_; }
  auto auto_range() const { return auto_range_; }
  auto min() const { return min_; }
  auto max() const { return max_; }
  const auto& title() const { return title_; }

 private:
  Float_t scale_;
  Float_t rel_error_;
  bool auto_range_;
  Float_t min_;
  Float_t max_;
  std::string title_;
};

class ScatterPlotOptions {
 public:
  ScatterPlotOptions(const JSONReader& params, const std::string& node);
  virtual ~ScatterPlotOptions() {};

  const auto& output_dir() const { return output_dir_; }

  const auto& draw_options() const { return draw_options_; }
  auto& title(std::string&& title) { title_ = title; return *this; }
  const auto& title() const { return title_; }

  auto& x() { return x_; }
  auto& y() { return y_; }
  const auto& x() const { return x_; }
  const auto& y() const { return y_; }

  const auto& marker_colors() const { return marker_colors_; }
  const auto& marker_sizes() const { return marker_sizes_; }
  const auto& marker_styles() const { return marker_styles_; }

 private:
  std::string output_dir_;

  std::string draw_options_;
  std::string title_;

  AxisOptions x_;
  AxisOptions y_;

  std::vector<Int_t> marker_colors_;
  std::vector<Float_t> marker_sizes_;
  std::vector<Int_t> marker_styles_;
};

#endif
