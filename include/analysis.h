#ifndef LUMIANALYSIS_INCLUDE_ANALYSIS_H_
#define LUMIANALYSIS_INCLUDE_ANALYSIS_H_

#include <string>
#include <vector>
#include <map>

#include "fit_results.h"
#include "single_run_data.h"

class Analysis {
 public:
  Analysis(std::string params_filepath);
  ~Analysis(){};

  int AnalyseTree(SingleRunData &this_run);
  int CreateAllRunPlots(const std::map<std::string, SingleRunData> &runs_data);
  int CreateLumiCurrentPlots(const SingleRunData &this_run);
  int CreateSingleRunPlots(const SingleRunData &this_run);
  int PrepareAnalysis(std::string params_filepath);
  int ReadCalibrations(std::string channels_filepath);
  int ReadChannelsList(std::string channels_list_filepath);
  void ReadParams(std::string params_filepath);
  int RunAnalysis();
  int WriteCurrentsToFile(std::string run_name);

  int CalcFCalLumi(SingleRunData &this_run);
  int CalcFCalMu(SingleRunData &this_run);
  int CreateBenedettoOutput(const SingleRunData &this_run) const;

 private:
  bool verbose_;

  bool do_Benedetto_;

  double f_rev_;
  double x_sec_;
  int ref_run_number_;
  double corr_A_;
  double corr_C_;
  double corr_Avg_;

  std::string params_filepath_;
  std::string calibrations_filepath_;
  std::string channels_list_filepath_;
  std::string pedestals_dir_;
  std::string trees_dir_;
  std::string run_list_dir_;
  std::string base_output_dir_;
  std::string plots_output_dir_;
  std::string fit_results_output_dir_;
  std::string calibrations_output_dir_;
  std::string geometric_analysis_output_dir_;
  std::string benedetto_output_dir_;

  std::vector<std::string> plot_types_;

  bool retrieve_timestamps_;
  bool retrieve_currents_;
  bool retrieve_lumi_BCM_;
  bool retrieve_lumi_FCal_;

  struct ChannelCalibration {
    Float_t slope;
    Float_t intercept;
  };

  std::map<std::string, ChannelCalibration> channels_list_;
};

#endif
