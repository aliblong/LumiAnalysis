#ifndef LUMIANALYSIS_INCLUDE_ANALYSIS_H_
#define LUMIANALYSIS_INCLUDE_ANALYSIS_H_

#include <string>
#include <vector>

using std::string;
using std::vector;

class Analysis {
 public:
  Analysis(string params_filepath);
  ~Analysis(){};
  int AnalyseTree(string tree_name);
  void ClearSingleRunVectors();
  int CreateSingleRunPlots(string run_name);
  int PrepareAnalysis(string params_filepath);
  int ReadCalibrations(string channels_filepath);
  int ReadChannelsList(string channels_list_filepath);
  int ReadPedestals(string pedestals_dir, string run_name);
  void ReadParams(string params_filepath);
  int RunAnalysis();
  int WriteCurrentsToFile(string run_name);

 private:
  bool verbose_;

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
  string output_dir_;

  vector<string> plot_types_;

  bool retrieve_currents_;
  bool retrieve_lumi_BCM_;
  bool retrieve_lumi_FCal_;

  struct Channel {
    Channel(string cn)
      : channel_name(cn),
        pedestal(0.),
        slope(0.),
        intercept(0.) {}
    string channel_name;
    float pedestal;
    float slope;
    float intercept;
  };

  vector<Channel> channels_list_;
  vector< vector<float> > currents_;
  vector<float> lumi_BCM_;
};

#endif
