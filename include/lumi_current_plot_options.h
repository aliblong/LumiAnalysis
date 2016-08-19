#ifndef LUMIANALYSIS_INCLUDE_LUMI_CURRENT_PLOTOPTIONS_H_
#define LUMIANALYSIS_INCLUDE_LUMI_CURRENT_PLOTOPTIONS_H_

#include <string>

#include "scatter_plot_options.h"

class LumiCurrentPlotOptions : public ScatterPlotOptions {
 public:
  LumiCurrentPlotOptions(const JSONReader& params, const std::string& node, bool verbose);
  ~LumiCurrentPlotOptions(){};

  auto& run_name(std::string&& run_name) { run_name_ = run_name; return *this; }
  auto& channel_name(std::string&& channel_name) { channel_name_ = channel_name; return *this; }

  const auto& run_name() const { return run_name_; }
  const auto& channel_name() const { return channel_name_; }

  const auto& plots_dir() const { return plots_dir_; }
  const auto& raw_fit_results_dir() const { return raw_fit_results_dir_; }
  const auto& calibration_results_dir() const { return calibration_results_dir_; }
  const auto& geometric_results_dir() const { return geometric_results_dir_; }

  auto do_individual() const { return do_individual_; }
  auto do_sum() const { return do_sum_; }
  auto do_fit() const { return do_fit_; }
  auto print_plots() const { return print_plots_; }

  auto fit_fix_intercept() const { return fit_fix_intercept_; }
  const auto& fit_options() const { return fit_options_; }
  auto fit_show_legend() const { return fit_show_legend_; }
  auto fit_line_color() const { return fit_line_color_; }
  auto fit_line_width() const { return fit_line_width_; }

 private:
  std::string run_name_;
  std::string channel_name_;

  std::string plots_dir_;
  std::string raw_fit_results_dir_;
  std::string calibration_results_dir_;
  std::string geometric_results_dir_;

  bool do_individual_;
  bool do_sum_;
  bool do_fit_;
  bool print_plots_;

  bool fit_fix_intercept_;
  std::string fit_options_;
  bool fit_show_legend_;
  int fit_line_color_;
  Float_t fit_line_width_;
};

#endif
