#ifndef LUMIANALYSIS_INCLUDE_FIT_RESULTS_H_
#define LUMIANALYSIS_INCLUDE_FIT_RESULTS_H_

#include <string>

#include <TF1.h>

#include <lumi_current_plot_options.h>

class FitResults {
 public:
  FitResults();
  FitResults(const TF1 &fit, double x_scale, double y_scale);
  ~FitResults();
  FitResults& operator=(const FitResults &rhs);

  void FromFit(const TF1 &fit, double x_scale, double y_scale);

  Float_t slope = 0.0;
  Float_t slope_err = 0.0;
  Float_t intercept = 0.0;
  Float_t intercept_err = 0.0;
  Float_t chi_squared = 0.0;
  UInt_t nDoF = 0;
  Bool_t is_short = false;

  Float_t calibration_slope = 0.0;
  Float_t calibration_intercept = 0.0;
};

#endif
