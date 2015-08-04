#ifndef LUMIANALYSIS_INCLUDE_MU_DEP_PLOTOPTIONS_H_
#define LUMIANALYSIS_INCLUDE_MU_DEP_PLOTOPTIONS_H_

#include <string>

#include "scatter_plot_options.h"

class MuDepPlotOptions : public ScatterPlotOptions {
 public:
  MuDepPlotOptions(const std::string& params_filepath);
  ~MuDepPlotOptions(){};

  const std::string& plots_dir() const { return plots_dir_; }

  const std::string& draw_options() const { return draw_options_; }

  auto& title(std::string&& title) { title_ = title; return *this; }
  const std::string& title() const { return title_; }

  int marker_color() const { return marker_color_; }
  Float_t marker_size() const { return marker_size_; }
  int marker_style() const { return marker_style_; }

  Float_t x_scale() const { return x_scale_; }
  Float_t x_rel_error() const { return x_rel_error_; }
  bool x_auto_range() const { return x_auto_range_; }
  Float_t x_min() const { return x_min_; }
  Float_t x_max() const { return x_max_; }
  const std::string& x_title() const { return x_title_; }

  Float_t y_scale() const { return y_scale_; }
  Float_t y_rel_error() const { return y_rel_error_; }
  bool y_auto_range() const { return y_auto_range_; }
  Float_t y_min() const { return y_min_; }
  Float_t y_max() const { return y_max_; }
  const std::string& y_title() const { return y_title_; }

 private:
  std::string plots_dir_;

  bool fit_fix_intercept_;
  std::string fit_options_;
  bool fit_show_legend_;
  int fit_line_color_;
  Float_t fit_line_width_;

  std::string draw_options_;
  std::string title_;

  int marker_color_;
  Float_t marker_size_;
  int marker_style_;

  Float_t x_scale_;
  Float_t x_rel_error_;
  bool x_auto_range_;
  Float_t x_min_;
  Float_t x_max_;
  std::string x_title_;

  Float_t y_scale_;
  Float_t y_rel_error_;
  bool y_auto_range_;
  Float_t y_min_;
  Float_t y_max_;
  std::string y_title_;
};

#endif
