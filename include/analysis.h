#ifndef LUMIANALYSIS_INCLUDE_ANALYSIS_H_
#define LUMIANALYSIS_INCLUDE_ANALYSIS_H_

#include <string>
#include <vector>

#include "boost/container/flat_map.hpp"
#include "boost/optional.hpp"

#include "detector.h"
#include "error.h"
#include "fit_results.h"
#include "void.h"

class Run;

// Contains analysis-wide parameters and methods which constitute the control
//   flow of the analysis
class Analysis {
 public:
  Analysis(std::string&& params_filepath);
  ~Analysis(){};

  struct ChannelCalibration {
    Float_t slope;
    Float_t intercept;
  };

  void CreateAllRunPlots(const boost::container::flat_map<std::string, Run> &runs_data);
  Error::Expected<Void> PrepareAnalysis();
  static Error::Expected<Void> ReadCalibrations(
      boost::container::flat_map<std::string, ChannelCalibration>* channel_calibrations,
      std::string primary_calibrations_filepath
  );
  Error::Expected<Void> ReadChannels();
  Error::Expected<Void> RunAnalysis();

  Detector::Name detector();
  bool verbose();

  double f_rev();
  double x_sec();

  int ref_run_number();
  double anchoring_factor_A();
  double anchoring_factor_C();
  double anchoring_factor_Avg();

  std::string primary_calibrations_filepath();
  std::string calibrations_dir();
  std::string channels_list_filepath();
  std::string pedestals_dir();
  std::string trees_dir();
  std::string currents_dir();
  std::string run_list_dir();

  std::vector<std::string> plot_types();

  std::string reference_lumi_algo();

  bool apply_LUCID_mu_corr();

  bool use_beamspot_corr();
  std::vector<Float_t> beamspot_corr_params();

  bool use_start_of_fill_pedestals();
  bool use_baseline_subtraction_from_fit();
  bool do_benedetto();
  std::string benedetto_output_dir();

  const auto& channel_calibrations() const { return channel_calibrations_; }

  const boost::container::flat_map<std::string, int>& n_bunches_from_file();
  const boost::container::flat_map<std::string, std::vector<int>>& custom_LB_bounds();

  const auto& params() const { return params_; }

  auto retrieve_timestamps() const { return retrieve_timestamps_; }
  auto retrieve_currents() const { return retrieve_currents_; }
  auto retrieve_lumi_ofl() const { return retrieve_lumi_ofl_; }
  auto retrieve_lumi_LAr() const { return retrieve_lumi_LAr_; }
  auto retrieve_mu_LAr() const { return retrieve_mu_LAr_; }
  auto retrieve_beamspot() const { return retrieve_beamspot_; }

 private:
  JSONReader params_;

  boost::optional<Detector::Name> detector_;
  boost::optional<bool> verbose_;

  boost::optional<double> f_rev_;
  boost::optional<double> x_sec_;

  boost::optional<int> ref_run_number_;
  boost::optional<double> anchoring_factor_A_;
  boost::optional<double> anchoring_factor_C_;
  boost::optional<double> anchoring_factor_Avg_;

  boost::optional<std::string> primary_calibrations_filepath_;
  boost::optional<std::string> calibrations_dir_;
  boost::optional<std::string> channels_list_filepath_;
  boost::optional<std::string> pedestals_dir_;
  boost::optional<std::string> trees_dir_;
  boost::optional<std::string> currents_dir_;
  boost::optional<std::string> run_list_dir_;

  boost::optional<std::vector<std::string>> plot_types_;

  boost::optional<std::string> reference_lumi_algo_;

  boost::optional<bool> apply_LUCID_mu_corr_;

  boost::optional<bool> use_beamspot_corr_;
  boost::optional<std::vector<Float_t>> beamspot_corr_params_;

  bool retrieve_timestamps_ = false;
  bool retrieve_currents_ = false;
  bool retrieve_lumi_ofl_ = false;
  bool retrieve_lumi_LAr_ = false;
  bool retrieve_mu_LAr_ = false;
  bool retrieve_beamspot_ = false;

  boost::optional<bool> use_start_of_fill_pedestals_;
  boost::optional<bool> use_baseline_subtraction_from_fit_;

  boost::optional<bool> do_benedetto_;
  boost::optional<std::string> benedetto_output_dir_;

  boost::container::flat_map<std::string, ChannelCalibration> channel_calibrations_;
  boost::optional<boost::container::flat_map<std::string, std::vector<int>>> custom_LB_bounds_;
  boost::optional<boost::container::flat_map<std::string, int>> n_bunches_from_file_;
};

#endif
