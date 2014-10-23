#ifndef LUMIANALYSIS_INCLUDE_LUMI_CURRENT_PLOTOPTIONS_H_
#define LUMIANALYSIS_INCLUDE_LUMI_CURRENT_PLOTOPTIONS_H_

#include <string>

class LumiCurrentPlotOptions {
 public:
  LumiCurrentPlotOptions(std::string params_filepath);
  ~LumiCurrentPlotOptions(){};

 //private:
  bool do_fit;

  std::string fit_options;
  bool fit_show_legend;
  int fit_line_color;
  float fit_line_width;

  std::string draw_options;

  int marker_color;
  float marker_size;
  int marker_style;

  float x_scale;
  bool x_auto_range;
  float x_min;
  float x_max;
  std::string x_title;

  float y_scale;
  bool y_auto_range;
  float y_min;
  float y_max;
  std::string y_title;
};

#endif
