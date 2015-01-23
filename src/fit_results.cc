#include <fit_results.h>

FitResults::FitResults()
  : slope(0.0),
    slope_err(0.0),
    intercept(0.0),
    intercept_err(0.0),
    chi_squared(0.0),
    nDoF(0),
    is_short(false),
    calibration_slope(0.0),
    calibration_intercept(0.0) {}

FitResults::FitResults(const TF1 &fit,
                       const LumiCurrentPlotOptions &plot_options)
  : FitResults() {

  this->FromFit(fit, plot_options);
}

FitResults::~FitResults() {}

FitResults& FitResults::operator=(const FitResults &rhs) {
  slope = rhs.slope;
  slope_err = rhs.slope_err;
  intercept = rhs.intercept;
  intercept_err = rhs.intercept_err;
  chi_squared = rhs.chi_squared;
  nDoF = rhs.nDoF;
  is_short = rhs.is_short;
  calibration_slope = rhs.calibration_slope;
  calibration_intercept = rhs.calibration_intercept;
  return *this;
}

void FitResults::FromFit(const TF1 &fit,
                         const LumiCurrentPlotOptions &plot_options) {
  slope = fit.GetParameter(1);
  slope_err = fit.GetParError(1);
  intercept = fit.GetParameter(0);
  intercept_err = fit.GetParError(0);
  chi_squared = fit.GetChisquare();
  nDoF = fit.GetNDF();
  // Convert to calibration used to get lumi = calib*current
  //   current*y_scale = m*lumi*x_scale + intercept
  //   lumi = 1 / (x_scale*m) * (current*y_scale + intercept)
  calibration_slope = plot_options.y_scale() / (slope*plot_options.x_scale());
  calibration_intercept = (-1) * intercept / (slope*plot_options.x_scale());
}
