#ifndef LUMIANALYSIS_INCLUDE_FIT_RESULTS_H_
#define LUMIANALYSIS_INCLUDE_FIT_RESULTS_H_

#include <string>

struct FitResults {
  /*
  void Set(string CN, float m, float m_E, float b,
           float b_E, float X2, int ndof) {
    channel_name = CN;
    slope = m;
    slope_error = m_E;
    intercept = b;
    intercept_error = b_E;
    chi_squared = X2;
    nDoF = ndof;
  }
  */
  std::string channel_name;
  float slope;
  float slope_err;
  float intercept;
  float intercept_err;
  float chi_squared;
  int nDoF;
};

#endif