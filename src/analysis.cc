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
#include "plot_options.h"

using std::string;
using std::vector;

using std::cout;
using std::cerr;
using std::endl;

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

void DumpVector(const vector<float> &vec) {
  for (auto element: vec) cout << element << endl;
}

bool IsZeroes(const vector<float> &vec) {
  float epsilon = 0.0000000001;
  for (auto element: vec) {
    if (element > epsilon) return false;
  }
  return true;
}

} // Unnamed namespace

Analysis::Analysis(string params_filepath) {
  params_filepath_ = params_filepath;
  int err = PrepareAnalysis(params_filepath);
  if (err) {
    cerr << "ERROR: in PrepareAnalysis(string params_filepath)" << endl;
  }
}

int Analysis::AnalyseTree(string run_name) {
  cout << "Analysing sample " << run_name << endl;
  ClearSingleRunVectors();
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

void Analysis::ClearSingleRunVectors() {
  currents_.clear();
  lumi_BCM_.clear();
}

int Analysis::CreateSingleRunPlots(string run_name) {
  cout << "Creating plots for sample " << run_name << endl;
  if (plot_types_.size() == 0) {
    cerr << "ERROR: no plots types selected to plot" << endl;
    return 1;
  }

  Plotter plotter;
  for (auto plot_type: plot_types_) {
    if (plot_type == "lumi_current") {
      cout << "    " << "Making lumi vs. current plots" << endl;
      PlotOptions plot_options(params_filepath_, plot_type);
      vector<FitResults> fit_results;
      int num_channels = used_channels_.size();
      for (int iChannel=0; iChannel < num_channels; ++iChannel) {
        auto this_channel_current = currents_.at(iChannel);
        auto this_channel_name = used_channels_.at(iChannel).channel_name;
        if ( IsZeroes(this_channel_current) ) {
          cout << "Skipping channel with zero current: " << this_channel_name
               << endl;
          continue;
        }
        FitResults this_channel_fit_results;
        int err_plot = plotter.PlotLumiCurrent(lumi_BCM_,
                                               this_channel_current,
                                               run_name,
                                               this_channel_name,
                                               output_dir_,
                                               plot_options,
                                               this_channel_fit_results);
        if (err_plot) {
          cerr << "ERROR: in Plotter::PlotLumiCurrent()" << endl;
          return 2;
        }
        fit_results.push_back(this_channel_fit_results);
      }
      int err_fit_results = plotter.SaveFitResults(fit_results,
                                                   run_name,
                                                   output_dir_);
    }
  }
  return 0;
}

int Analysis::PrepareAnalysis(string params_filepath) {
  ReadParams(params_filepath);
  for (auto plot_type: plot_types_) {
    if (plot_type == "lumi_current" ) {
      retrieve_currents_ = true;
      retrieve_lumi_BCM_ = true;
    }
  }
  int err_RCF = ReadCalibrations(calibrations_filepath_);
  if (err_RCF) {
    cerr << "ERROR: in ReadCalibrations(string calibrations_filepath"
         << ", string pedestals_dir)" << endl;
    return 1;
  }

  return 0;
}

int Analysis::ReadCalibrations(string calibrations_filepath) {
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

  return 0;
}

int Analysis::ReadPedestals(string pedestals_dir, string run_name) {
  for (auto &channel: used_channels_) {
    auto pedestals_filepath = pedestals_dir + run_name + ".dat";
    std::ifstream pedestals_file(pedestals_filepath);
    if (!pedestals_file) {
      cerr << "ERROR: Could not locate pedestals file \'"
           << pedestals_filepath << "\'" << endl;
      return 1;
    }

    string channel_name;
    float pedestal;
    float something;
    bool found_channel = false;
    while (pedestals_file >> channel_name >> pedestal >> something) {
      if (channel.channel_name == channel_name) {
        channel.pedestal = pedestal;
        found_channel = true;
        break;
      }
    }
    pedestals_file.close();

    if (!found_channel) {
      cerr << "ERROR: Could not locate pedestal for run " << run_name
           << ", channel " << channel.channel_name << endl;
      return 2;
    }
    found_channel = false;
  } // Used channels loop
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
  pedestals_dir_ = parameter_file.get<string>("filepaths.pedestals");
  trees_dir_ = parameter_file.get<string>("filepaths.trees");
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
  while ( getline(run_names_file, run_name) ) {
    if (run_name[0] == '#') {
      cout << "Skipping run " << run_name << endl;
      continue;
    }
    int err_pedestal = ReadPedestals(pedestals_dir_, run_name);
    if (err_pedestal) {
      cerr << "ERROR: In ReadPedestals(string pedestals_dir, string run_name)"
           << endl;
      return 1;
    }
    int err_analyse = AnalyseTree(run_name);
    if (err_analyse) {
      cerr << "ERROR: In AnalyseTree(string run_name)" << endl;
      return 1;
    }
    int err_plots = CreateSingleRunPlots(run_name);
    if (err_plots) {
      cerr << "ERROR: In CreateSingleRunPlots(string run_name)" << endl;
      return 1;
    }
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
