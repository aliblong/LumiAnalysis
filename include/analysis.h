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
  void ClearVectors();
  int CreateSingleRunPlots(string run_name);
  int PrepareAnalysis(string params_filepath);
  int ReadChannelsCalibAndPed(string channels_filepath,
                              string pedestals_filepath);
  void ReadParams(string params_filepath);
  int RunAnalysis();
  int WriteCurrentsToFile(string run_name);

 private:
  double f_rev_;
  double x_sec_;
  int ref_run_number_;
  double corr_A_;
  double corr_C_;
  double corr_Avg_;
  vector<string> plot_types_;
  string params_filepath_;
  string calibrations_filepath_;
  string pedestals_filepath_;
  string trees_dir_;
  string baselines_dir_;
  string run_list_dir_;
  string output_dir_;

  bool retrieve_currents_;
  bool retrieve_lumi_BCM_;
  struct ChannelCalibration {
    ChannelCalibration(string mn, float s, float i) {
      channel_name = mn;
      slope = s;
      intercept = i;
    }
    float pedestal;
    string channel_name;
    float slope;
    float intercept;
  };

  vector<ChannelCalibration> used_channels_;
  vector< vector<float> > currents_;
  vector<float> lumi_BCM_;
};

#endif
