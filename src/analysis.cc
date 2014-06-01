#include <iostream>
#include <string>
#include <vector>
#include <fstream>

#include "TFile.h"
#include "TTree.h"

#include "analysis.h"
#include "branch_array_buffer_sizes.h"
#include "json_reader.h"
#include "plotter.h"

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::vector;

namespace {

template <typename T>
int WriteVectorToFile(const vector<T> &vec, string filepath) {
  std::ofstream out_file(filepath);
  if (!out_file.is_open()) {
    cerr << "ERROR: could not open file \'" << filepath << "\'" << endl;
    return 1;
  }

  for (auto iElement: vec) out_file << iElement << '\n';
  out_file.close();
  return 0;
}

int GetNumLumiBlocks(TTree *tree) {
  int num_lumi_blocks = 0;
  tree->SetBranchAddress("NLB", &num_lumi_blocks);
  tree->GetEntry(0);
  return num_lumi_blocks;
}

} // Unnamed namespace

Analysis::Analysis(string params_filepath) {
  int err = PrepareAnalysis(params_filepath);
  if (err) {
    cerr << "ERROR: in PrepareAnalysis(string params_filepath)" << endl;
  }
}

int Analysis::AnalyseTree(string run_name) {
  string filepath = trees_dir_+run_name+".root";
  TFile *this_file = TFile::Open(filepath.c_str());
  if (!this_file) {
    cerr << "ERROR: Could not locate file: " << filepath << endl;
    return 1;
  }

  string tree_name = "t_"+run_name;
  TTree *this_tree = static_cast<TTree*>( this_file->Get(tree_name.c_str()) );
  if (!this_tree) {
    cerr << "ERROR: tree \'" << tree_name << "\' not found in file \'" <<
            filepath << "\'" << endl;
    return 2;
  }

  int num_lumi_blocks = GetNumLumiBlocks(this_tree);
  int num_channels = used_channels_.size();

  float currents[gMaxNumChannels][gMaxNumLB];
  float lumi_BCM[gMaxNumLB];
  int quality[gMaxNumLB];
  std::unique_ptr< vector<string> > beam_mode(new vector<string>);

  if (retrieve_currents_) {
    for (int iChannel=0; iChannel < num_channels; ++iChannel) {
      string branch_name = "current_"+used_channels_[iChannel].channel_name;
      this_tree->SetBranchAddress((branch_name).c_str(), &currents[iChannel]);
    }
  }
  if (retrieve_lumi_BCM_) {
    this_tree->SetBranchAddress("ofl_lumi_pref",&lumi_BCM);
  }

  this_tree->SetBranchAddress("quality", &quality);
  this_tree->SetBranchAddress("mode", &beam_mode);

  this_tree->GetEntry(0);

  if (retrieve_currents_) {
    for (int iChannel=0; iChannel < num_channels; ++iChannel) {
      currents_.emplace_back();
      float this_channel_pedestal = used_channels_.at(iChannel).pedestal;

      for (int iLumiBlock=0; iLumiBlock < num_lumi_blocks; ++iLumiBlock) {
        if (quality[iLumiBlock] && beam_mode->at(iLumiBlock) == "STABLE BEAMS") {
          currents_.at(iChannel).push_back( currents[iChannel][iLumiBlock] -
                                            this_channel_pedestal );
        }
      }
    }
  }

  if (retrieve_lumi_BCM_) {
    for (int iLumiBlock=0; iLumiBlock < num_lumi_blocks; ++iLumiBlock) {
      if (quality[iLumiBlock] && beam_mode->at(iLumiBlock) == "STABLE BEAMS") {
        lumi_BCM_.push_back( lumi_BCM[iLumiBlock] );
      }
    }
  }

  return 0;
}

int Analysis::CreateSingleRunPlots(string run_name) {
  if (plot_types_.size() == 0) {
    cerr << "ERROR: no plots types selected to plot" << endl;
    return 1;
  }

  Plotter plotter;
  for (auto plot_type: plot_types_) {
    if (plot_type == "lc") {
      int num_channels = used_channels_.size();
      for (int iChannel=0; iChannel < num_channels; ++iChannel) {
        auto *this_channel_current = &currents_.at(iChannel);
        auto this_channel_name = used_channels_.at(iChannel).channel_name;
        int err_plot = plotter.PlotLumiCurrent(lumi_BCM_,
                                               *this_channel_current,
                                               run_name,
                                               this_channel_name,
                                               output_dir_);
        if (err_plot) {
          cerr << "ERROR: in Plotter::PlotLumiCurrent()" << endl;
          return 2;
        }
      }
    }
  }
  return 0;
}

int Analysis::PrepareAnalysis(string params_filepath) {
  ReadParams(params_filepath);
  for (auto plot_type: plot_types_) {
    if (plot_type == "lc" ) {
      retrieve_currents_ = true;
      retrieve_lumi_BCM_ = true;
    }
  }
  int err_RCF = ReadChannelsCalibAndPed(calibrations_filepath_,
                                        pedestals_filepath_);
  if (err_RCF) {
    cerr << "ERROR: in ReadChannelsCalibAndPed(string calibrations_filepath"
         << ", string pedestals_filepath)" << endl;
    return 1;
  }

  return 0;
}

int Analysis::ReadChannelsCalibAndPed(string calibrations_filepath,
                                      string pedestals_filepath) {
  std::ifstream calibrations_file(calibrations_filepath);
  if (!calibrations_file) {
    cerr << "ERROR: Could not locate calibrations file \'"
         << calibrations_filepath << "\'" << endl;
    return 1;
  }

  string channel_name;
  float slope;
  float intercept;
  while (calibrations_file >> channel_name >> slope >> intercept) {
    used_channels_.emplace_back(channel_name, slope, intercept);
  }
  calibrations_file.close();

  std::ifstream pedestals_file(pedestals_filepath);
  if (!pedestals_file) {
    cerr << "ERROR: Could not locate pedestals file \'"
         << calibrations_filepath << "\'" << endl;
    return 1;
  }

  float pedestal;
  while (pedestals_file >> channel_name >> pedestal) {
    for (auto &channel: used_channels_) {
      if (channel.channel_name == channel_name) {
        channel.pedestal = pedestal;
        break;
      }
    }
  }
  pedestals_file.close();

  return 0;
}

void Analysis::ReadParams(string params_filepath) {
  JSONReader parameter_file(params_filepath);

  // Used to calculate
  f_rev_ = parameter_file.get<double>("lumi_calculation.f_rev");
  x_sec_ = parameter_file.get<double>("lumi_calculation.x_sec");

  // Value of 0 for reference run means set no corrections
  ref_run_number_ = parameter_file.get<int>("ref_run_number");
  string ref_run_str = std::to_string(ref_run_number_);
  corr_A_ = parameter_file.get<double>("corrections."+ref_run_str+".A");
  corr_C_ = parameter_file.get<double>("corrections."+ref_run_str+".C");
  corr_Avg_ = parameter_file.get<double>("corrections."+ref_run_str+".Avg");

  // plot_types_ was initialized in Analysis ctor, so can't use c++11 move
  //   semantics; must instead use swaptimization
  parameter_file.get_vector<string>("plot_types").swap(plot_types_);

  // Relevant directories
  calibrations_filepath_ = parameter_file.get<string>("filepaths.calibrations");
  pedestals_filepath_ = parameter_file.get<string>("filepaths.pedestals");
  trees_dir_ = parameter_file.get<string>("filepaths.trees");
  baselines_dir_ = parameter_file.get<string>("filepaths.baselines");
  run_list_dir_ = parameter_file.get<string>("filepaths.run_list");
  output_dir_ = parameter_file.get<string>("filepaths.output");
}

int Analysis::RunAnalysis() {
  ifstream run_names_file(run_list_dir_);
  if (!run_names_file) {
    cerr << "ERROR: Could not locate file: " << run_list_dir_ << endl;
    return 1;
  }

  string run_name;
  while (run_names_file >> run_name) {
    int err_analyse = AnalyseTree(run_name);
    if (err_analyse) {
      cerr << "ERROR: In AnalyseTree(string run_name)" << endl;
      return 1;
    }
    int err_plots = CreateSingleRunPlots(run_name);
    //int err_write = WriteVectorToFile(lumi_BCM_, "testout.dat");
    /*
    string out_filepath_lumi_BCM = out_dir_+"lumi_BCM/"+run_name+".dat";
    int err_lumi_BCM_write = WriteVectorToFile(lumi_BCM,
                                               out_filepath_lumi_BCM);
    if (err_lumi_BCM_write) {
      cerr << "ERROR: In WriteVectorToFile(lumi_BCM_, "
           << run_name << ")" << endl;
      return 1;
    }
    */
  }
  return 0;
}

/*
int WriteCurrentsToFile(vector<float> currents, string run_name) {
  std::ofstream out_file(out_filepath);
  if (!out_file.is_open()) {
    cerr << "ERROR: Could not open file \'" << out_filepath << "\'" << endl;
    return 1;
  }
  for (auto iEvent: currents_) {
    for (auto iCurrent: iEvent) out_file << iCurrent << ' ';
    out_file << '\n';
  }
  out_file.close();
  return 0;
}
*/
