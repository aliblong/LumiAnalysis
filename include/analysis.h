#ifndef LUMIANALYSIS_INCLUDE_ANALYSIS_H_
#define LUMIANALYSIS_INCLUDE_ANALYSIS_H_

#include <string>
#include <vector>

#include "boost/container/flat_map.hpp"

#include "detector.h"
#include "error.h"
#include "fit_results.h"
#include "void.h"

class SingleRunData;

// Contains analysis-wide parameters and methods which constitute the control
//   flow of the analysis
class Analysis {
 public:
  Analysis(std::string&& params_filepath);
  ~Analysis(){};

  void CreateAllRunPlots(const boost::container::flat_map<std::string, SingleRunData> &runs_data);
  Error::Expected<Void> PrepareAnalysis();
  Error::Expected<Void> ReadCalibrations();
  Error::Expected<Void> ReadChannels();
  Error::Expected<Void> ReadParams();
  Error::Expected<Void> RunAnalysis();

  auto detector() const { return detector_; }
  auto verbose() const { return verbose_; }

  auto f_rev() const { return f_rev_; }
  auto x_sec() const { return x_sec_; }
  auto ref_run_number() const { return ref_run_number_; }
  auto corr_A() const { return corr_A_; }
  auto corr_C() const { return corr_C_; }
  auto corr_Avg() const { return corr_Avg_; }

  const auto& params_filepath() const { return params_filepath_; }
  const auto& calibrations_filepath() const { return calibrations_filepath_; }
  const auto& channels_list_filepath() const { return channels_list_filepath_; }
  const auto& pedestals_dir() const { return pedestals_dir_; }
  const auto& trees_dir() const { return trees_dir_; }
  const auto& run_list_dir() const { return run_list_dir_; }

  const auto& plot_types() const { return plot_types_; }

  auto retrieve_timestamps() const { return retrieve_timestamps_; }
  auto retrieve_currents() const { return retrieve_currents_; }
  auto retrieve_lumi_ofl() const { return retrieve_lumi_ofl_; }
  auto retrieve_lumi_FCal() const { return retrieve_lumi_FCal_; }

  auto use_start_of_fill_pedestals() const { return use_start_of_fill_pedestals_; }

  auto do_benedetto() const { return do_benedetto_; }
  const auto& benedetto_output_dir() const { return benedetto_output_dir_; }

  const auto& channel_calibrations() const { return channel_calibrations_; }
  const auto& custom_LB_bounds() const { return custom_LB_bounds_; }

 private:
  Detector::Name detector_ = Detector::Name::FCal;
  bool verbose_ = false;

  double f_rev_ = 0.0;
  double x_sec_ = 0.0;

  int ref_run_number_ = 0;
  double corr_A_ = 0.0;
  double corr_C_ = 0.0;
  double corr_Avg_ = 0.0;

  std::string params_filepath_;
  std::string calibrations_filepath_;
  std::string channels_list_filepath_;
  std::string pedestals_dir_;
  std::string trees_dir_;
  std::string run_list_dir_;

  std::vector<std::string> plot_types_;

  bool retrieve_timestamps_ = false;
  bool retrieve_currents_ = false;
  bool retrieve_lumi_ofl_ = false;
  bool retrieve_lumi_FCal_ = false;

  bool use_start_of_fill_pedestals_ = false;

  bool do_benedetto_ = false;
  std::string benedetto_output_dir_;

  struct ChannelCalibration {
    Float_t slope;
    Float_t intercept;
  };

  boost::container::flat_map<std::string, ChannelCalibration> channel_calibrations_;
  boost::container::flat_map<std::string, std::vector<int>> custom_LB_bounds_;
};

#endif
