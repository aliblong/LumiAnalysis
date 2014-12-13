#ifndef LUMIANALYSIS_INCLUDE_FIT_RESULTS_H_
#define LUMIANALYSIS_INCLUDE_FIT_RESULTS_H_

#include <string>

#include <TF1.h>

#include <lumi_current_plot_options.h>

class FitResults {
 public:
  FitResults();
  ~FitResults();
  FitResults& operator=(const FitResults &rhs);

  void FromFit(const TF1 *fit, const LumiCurrentPlotOptions &plot_options);

  Float_t slope;
  Float_t slope_err;
  Float_t intercept;
  Float_t intercept_err;
  Float_t chi_squared;
  UInt_t nDoF;
  Bool_t is_short;

  Float_t calibration_slope;
  Float_t calibration_intercept;
};

#endif
