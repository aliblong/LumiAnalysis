#ifndef LUMIANALYSIS_INCLUDE_MU_DEP_PLOTOPTIONS_H_
#define LUMIANALYSIS_INCLUDE_MU_DEP_PLOTOPTIONS_H_

#include <string>

class MuDepPlotOptions {
 public:
  MuDepPlotOptions(const std::string& params_filepath);
  ~MuDepPlotOptions(){};

  const auto& plots_dir() const { return plots_dir_; }

  const auto& draw_options() const { return draw_options_; }

  auto marker_color() const { return marker_color_; }
  auto marker_size() const { return marker_size_; }
  auto marker_style() const { return marker_style_; }

  auto x_scale() const { return x_scale_; }
  auto x_rel_error() const { return x_rel_error_; }
  auto x_auto_range() const { return x_auto_range_; }
  auto x_min() const { return x_min_; }
  auto x_max() const { return x_max_; }
  const auto& x_title() const { return x_title_; }

  auto y_scale() const { return y_scale_; }
  auto y_rel_error() const { return y_rel_error_; }
  auto y_auto_range() const { return y_auto_range_; }
  auto y_min() const { return y_min_; }
  auto y_max() const { return y_max_; }
  const auto& y_title() const { return y_title_; }

 private:
  std::string plots_dir_;

  bool fit_fix_intercept_;
  std::string fit_options_;
  bool fit_show_legend_;
  int fit_line_color_;
  Float_t fit_line_width_;

  std::string draw_options_;

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
