#ifndef LUMIANALYSIS_INCLUDE_PLOTOPTIONS_H_
#define LUMIANALYSIS_INCLUDE_PLOTOPTIONS_H_

#include <string>

class PlotOptions {
  typedef typename std::string string;
 public:
  PlotOptions(string params_filepath, string plot_type);
  ~PlotOptions(){};

// private:
  bool do_fit;

  string fit_options;
  bool fit_show_legend;
  int fit_line_color;
  float fit_line_width;

  string draw_options;

  int marker_color;
  int marker_style;
  float marker_size;

  float x_scale;
  bool x_auto_range;
  float x_min;
  float x_max;
  string x_title;

  float y_scale;
  bool y_auto_range;
  float y_min;
  float y_max;
  string y_title;
};

#endif
