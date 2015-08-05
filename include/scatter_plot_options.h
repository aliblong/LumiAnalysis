#ifndef LUMIANALYSIS_INCLUDE_SCATTER_PLOT_PLOTOPTIONS_H_
#define LUMIANALYSIS_INCLUDE_SCATTER_PLOT_PLOTOPTIONS_H_

#include <string>

#include "Rtypes.h"

class ScatterPlotOptions {
 public:
  virtual ~ScatterPlotOptions() {};

  virtual const std::string& draw_options() const = 0;

  virtual const std::string& title() const = 0;

  virtual int marker_color() const = 0;
  virtual Float_t marker_size() const = 0;
  virtual int marker_style() const = 0;

  virtual Float_t x_scale() const = 0;
  virtual Float_t x_rel_error() const = 0;
  virtual bool x_auto_range() const = 0;
  virtual Float_t x_min() const = 0;
  virtual Float_t x_max() const = 0;
  virtual const std::string& x_title() const = 0;

  virtual Float_t y_scale() const = 0;
  virtual Float_t y_rel_error() const = 0;
  virtual bool y_auto_range() const = 0;
  virtual Float_t y_min() const = 0;
  virtual Float_t y_max() const = 0;
  virtual const std::string& y_title() const = 0;
};

#endif
