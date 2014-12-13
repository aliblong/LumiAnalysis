#ifndef LUMIANALYSIS_INCLUDE_MU_STAB_PLOTOPTIONS_H_
#define LUMIANALYSIS_INCLUDE_MU_STAB_PLOTOPTIONS_H_

#include <string>

#include "Rtypes.h"

class MuStabPlotOptions {
 public:
  MuStabPlotOptions(std::string params_filepath);
  ~MuStabPlotOptions(){};

  std::string draw_options() const { return draw_options_; }

  int marker_color_A() const { return marker_color_A_; }
  Float_t marker_size_A() const { return marker_size_A_; }
  int marker_style_A() const { return marker_style_A_; }
  int marker_color_C() const { return marker_color_C_; }
  Float_t marker_size_C() const { return marker_size_C_; }
  int marker_style_C() const { return marker_style_C_; }
  int marker_color_Avg() const { return marker_color_Avg_; }
  Float_t marker_size_Avg() const { return marker_size_Avg_; }
  int marker_style_Avg() const { return marker_style_Avg_; }

  bool x_auto_range() const { return x_auto_range_; }
  UInt_t low_bin() const { return low_bin_; }
  UInt_t high_bin() const { return high_bin_; }
  UInt_t n_bins() const { return n_bins_; }
  std::string x_title() const { return x_title_; }

  bool y_auto_range() const { return y_auto_range_; }
  Float_t y_min() const { return y_min_; }
  Float_t y_max() const { return y_max_; }
  std::string y_title() const { return y_title_; }

 private:
  std::string draw_options_;

  int marker_color_A_;
  Float_t marker_size_A_;
  int marker_style_A_;
  int marker_color_C_;
  Float_t marker_size_C_;
  int marker_style_C_;
  int marker_color_Avg_;
  Float_t marker_size_Avg_;
  int marker_style_Avg_;

  bool x_auto_range_;
  UInt_t low_bin_;
  UInt_t high_bin_;
  UInt_t n_bins_;
  std::string x_title_;

  bool y_auto_range_;
  Float_t y_min_;
  Float_t y_max_;
  std::string y_title_;
};

#endif
