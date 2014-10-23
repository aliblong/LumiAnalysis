#ifndef LUMIANALYSIS_INCLUDE_MU_STAB_PLOTOPTIONS_H_
#define LUMIANALYSIS_INCLUDE_MU_STAB_PLOTOPTIONS_H_

#include <string>

#include "Rtypes.h"

class MuStabPlotOptions {
 public:
  MuStabPlotOptions(std::string params_filepath);
  ~MuStabPlotOptions(){};

 //private:
  std::string draw_options;

  int marker_color_A;
  float marker_size_A;
  int marker_style_A;
  int marker_color_C;
  float marker_size_C;
  int marker_style_C;
  int marker_color_Avg;
  float marker_size_Avg;
  int marker_style_Avg;

  bool x_auto_range;
  UInt_t low_bin;
  UInt_t high_bin;
  UInt_t bin_size;
  std::string x_title;

  bool y_auto_range;
  float y_min;
  float y_max;
  std::string y_title;
};

#endif
