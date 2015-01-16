#ifndef LUMIANALYSIS_INCLUDE_LUMI_CURRENT_PLOTOPTIONS_H_
#define LUMIANALYSIS_INCLUDE_LUMI_CURRENT_PLOTOPTIONS_H_

#include <string>

class LumiCurrentPlotOptions {
 public:
  LumiCurrentPlotOptions(std::string params_filepath);
  ~LumiCurrentPlotOptions(){};

  bool do_individual() const { return do_individual_; }
  bool do_sum() const { return do_sum_; }
  bool do_fit() const { return do_fit_; }

  bool fit_fix_intercept() const { return fit_fix_intercept_; }
  std::string fit_options() const { return fit_options_; }
  bool fit_show_legend() const { return fit_show_legend_; }
  int fit_line_color() const { return fit_line_color_; }
  Float_t fit_line_width() const { return fit_line_width_; }

  std::string draw_options() const { return draw_options_; }

  int marker_color() const { return marker_color_; }
  Float_t marker_size() const { return marker_size_; }
  int marker_style() const { return marker_style_; }

  Float_t x_scale() const { return x_scale_; }
  Float_t x_rel_error() const { return x_rel_error_; }
  bool x_auto_range() const { return x_auto_range_; }
  Float_t x_min() const { return x_min_; }
  Float_t x_max() const { return x_max_; }
  std::string x_title() const { return x_title_; }

  Float_t y_scale() const { return y_scale_; }
  Float_t y_rel_error() const { return y_rel_error_; }
  bool y_auto_range() const { return y_auto_range_; }
  Float_t y_min() const { return y_min_; }
  Float_t y_max() const { return y_max_; }
  std::string y_title() const { return y_title_; }

 private:
  bool do_individual_;
  bool do_sum_;
  bool do_fit_;

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
