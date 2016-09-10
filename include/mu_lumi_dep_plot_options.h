#ifndef LUMIANALYSIS_INCLUDE_MU_LUMI_DEP_PLOTOPTIONS_H_
#define LUMIANALYSIS_INCLUDE_MU_LUMI_DEP_PLOTOPTIONS_H_

#include <string>

#include "scatter_plot_options.h"

class MuLumiDepPlotOptions : public ScatterPlotOptions {
 public:
  MuLumiDepPlotOptions(const JSONReader& params, const std::string& node);
  ~MuLumiDepPlotOptions(){};

  const std::string& file_name() const { return file_name_; }
  //auto& file_name(std::string&& file_name) { file_name_ = file_name; return *this; }

  auto do_fit() const { return do_fit_; }
  auto fit_fix_intercept() const { return fit_fix_intercept_; }
  const auto& fit_options() const { return fit_options_; }
  auto fit_show_legend() const { return fit_show_legend_; }
  auto fit_line_color() const { return fit_line_color_; }
  auto fit_line_width() const { return fit_line_width_; }

 private:
  std::string file_name_;

  bool do_fit_;
  bool fit_fix_intercept_;
  std::string fit_options_;
  bool fit_show_legend_;
  int fit_line_color_;
  Float_t fit_line_width_;
};

#endif
