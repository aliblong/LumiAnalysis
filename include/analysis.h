#ifndef LUMIANALYSIS_INCLUDE_ANALYSIS_H_
#define LUMIANALYSIS_INCLUDE_ANALYSIS_H_

#include <string>
#include <vector>
#include <map>

#include "single_run_data.h"

class Analysis {
  typedef typename std::string string;

 public:
  Analysis(string params_filepath);
  ~Analysis(){};

  int AnalyseTree(SingleRunData &this_run);
  int CreateSingleRunPlots(const SingleRunData &this_run);
  int PrepareAnalysis(string params_filepath);
  int ReadCalibrations(string channels_filepath);
  int ReadChannelsList(string channels_list_filepath);
  void ReadParams(string params_filepath);
  int RunAnalysis();
  int WriteCurrentsToFile(string run_name);

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

  string params_filepath_;
  string calibrations_filepath_;
  string channels_list_filepath_;
  string pedestals_dir_;
  string trees_dir_;
  string run_list_dir_;
  string base_output_dir_;
  string plots_output_dir_;
  string benedetto_output_dir_;

  std::vector<string> plot_types_;

  bool retrieve_currents_;
  bool retrieve_lumi_BCM_;
  bool retrieve_lumi_FCal_;

  struct ChannelCalibration {
    float slope;
    float intercept;
  };

  std::map<string, ChannelCalibration> channels_list_;
};

#endif
