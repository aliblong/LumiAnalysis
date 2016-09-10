#include "mu_lumi_dep_plot_options.h"

#include <string>

#include "Rtypes.h"

#include "json_reader.h"

using std::string;

// Decline to handle errors from missing node in both input and default parameter trees
MuLumiDepPlotOptions::MuLumiDepPlotOptions(const JSONReader& params, const std::string& node) :
  ScatterPlotOptions(params, node)
{
  // Fit options
  do_fit_ = *params.get<bool>(node+"fit.do");
  string fit_verbosity_option;
  int fit_verbosity = *params.get<int>(node+"fit.verbose");
  if (fit_verbosity == 2) {
    fit_verbosity_option = "V";
  }
  else if (fit_verbosity == 1) {
    fit_verbosity_option = "";
  }
  else { // verbosity 0
    fit_verbosity_option = "Q";
  }
  string fit_show_option;
  if (*params.get<bool>(node+"fit.show")) {
    fit_show_option = "";
  }
  else {
    fit_show_option = "0";
  }
  string additional_options = *params.get<string>(node+"fit.options");
  fit_options_ = fit_verbosity_option + fit_show_option + additional_options;

  fit_fix_intercept_ = *params.get<bool>(node+"fit.fix_intercept");
  fit_show_legend_ = *params.get<bool>(node+"fit.show_legend");
  fit_line_color_ = *params.get<int>(node+"fit.line_color");
  fit_line_width_ = *params.get<Float_t>(node+"fit.line_width");

  file_name_ = *params.get<string>(node+"file_name");
}
