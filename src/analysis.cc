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
  // Writes each element of a vector to a file at the specified path;
  //   elements are separated by newlines.

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
  // Retrieves the quantity stored in the NLB variable in a TTree.

  int num_lumi_blocks = 0;
  tree->SetBranchAddress("NLB", &num_lumi_blocks);
  tree->GetEntry(0);

  return num_lumi_blocks;
}

void DumpVector(const vector<float> &vec) {
  // Prints newline-separated vector elements to cout.
  for (auto element: vec) cout << element << endl;
}

bool IsZeroes(const vector<float> &vec) {
  // Checks whether or not the vector has all elements less than epsilon

  float epsilon = 0.0000000001;
  for (auto element: vec) {
    if (element > epsilon) return false;
  }
  return true;
}

} // Anonymous namespace

Analysis::Analysis(string params_filepath) {
  // Initializes various quantities required to run the analysis

  params_filepath_ = params_filepath;
  int err = PrepareAnalysis(params_filepath);
  if (err) {
    cerr << "ERROR: in PrepareAnalysis(string params_filepath)" << endl;
  }
}

int Analysis::AnalyseTree(string run_name) {
  // Runs the analysis for a single sample.
  // Saves relevant results to member variables of the analysis object.

  if (verbose_) {
    cout << "Analysing sample " << run_name << endl;
  }

  // Clears data used in plots of single samples, which should have been used
  //   already at this point
  ClearSingleRunVectors();

  // Get tree object from the sample of choice
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

  // The number of luminosity readings in the run
  int num_lumi_blocks = GetNumLumiBlocks(this_tree);

  // Allocate buffers big enough to hold the data (can't use dynamic memory
  //   because of how TTrees work)
  float currents[gMaxNumChannels][gMaxNumLB];
  float lumi_BCM[gMaxNumLB];
  int quality[gMaxNumLB];
  std::unique_ptr< vector<string> > beam_mode(new vector<string>);

  // Set variables to retrieve based on which plot types are being produced
  if (retrieve_currents_) {
    //for (int iChannel=0; iChannel < num_channels; ++iChannel) {
    unsigned iChannel = 0;
    for (const auto &channel: channels_list_) {
      string branch_name = "current_"+channel.first;
      this_tree->SetBranchAddress((branch_name).c_str(), &currents[iChannel]);
      ++iChannel;
    }
  }
  if (retrieve_lumi_BCM_) {
    this_tree->SetBranchAddress("ofl_lumi_pref",&lumi_BCM);
  }

  this_tree->SetBranchAddress("quality", &quality);
  this_tree->SetBranchAddress("mode", &beam_mode);

  // Populate those variables which have been set to branches
  this_tree->GetEntry(0);

  // Save results to member variables
  if (retrieve_currents_) {
    //for (int iChannel=0; iChannel < num_channels; ++iChannel) {
    unsigned iChannel = 0;
    for (const auto &channel: channels_list_) {
      string this_channel_name = channel.first;

      vector<float> current_vec;
      currents_.insert(std::make_pair(this_channel_name, current_vec));

      float this_channel_pedestal = channel.second.pedestal;
      auto &this_channel_currents = currents_.at(this_channel_name);

      // If beam quality is good, save pedestal-subtracted current
      for (int iLumiBlock=0; iLumiBlock < num_lumi_blocks; ++iLumiBlock) {
        if (quality[iLumiBlock] && beam_mode->at(iLumiBlock) == "STABLE BEAMS") {
          auto this_current = currents[iChannel][iLumiBlock] -
                              this_channel_pedestal;
          this_channel_currents.push_back(this_current);
        }
      }

      ++iChannel;
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
// Create those plots which use data from only a single sample

  if (verbose_) cout << "Creating plots for sample " << run_name << endl;
  if (plot_types_.size() == 0) {
    cerr << "ERROR: no plots types selected to plot" << endl;
    return 1;
  }

  Plotter plotter;
  for (auto plot_type: plot_types_) {
    if (plot_type == "lumi_current") {
      if (verbose_) cout << "    " << "Making lumi vs. current plots" << endl;

      PlotOptions plot_options(params_filepath_, plot_type);
      std::map<string, FitResults> fit_results;
      //int num_channels = channels_list_.size();

      for (const auto &channel: channels_list_) {
        auto this_channel_name = channel.first;
        auto this_channel_current = currents_.at(this_channel_name);

        if (IsZeroes(this_channel_current) && verbose_) {
          cout << "Skipping channel with zero current: " << this_channel_name
               << endl;
          continue;
        }

        FitResults this_channel_fit_results;
        if (channels_list_.at(this_channel_name).pedestal > 20) {
          this_channel_fit_results.is_short = true;
        } else {
          this_channel_fit_results.is_short = false;
        }

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

        if (plot_options.do_fit) {
          fit_results.insert(std::make_pair(this_channel_name,
                                            this_channel_fit_results));
        }
      }

      if (plot_options.do_fit) {
        int err_fit_results = plotter.SaveFitResults(fit_results,
                                                     run_name,
                                                     output_dir_);
      }
    }
  }

  return 0;
}

int Analysis::PrepareAnalysis(string params_filepath) {
// Reads in the analysis parameters (json) and channel information (text). 
// Sets flags for which data to retrieve based on which plots to produce.

  ReadParams(params_filepath);
  for (auto plot_type: plot_types_) {
    if (plot_type == "lumi_current" ) {
      retrieve_currents_ = true;
      retrieve_lumi_BCM_ = true;
    }
  }

  int err_RCL = ReadChannelsList(channels_list_filepath_);
  if (err_RCL) {
    cerr << "ERROR: in ReadCalibrations(string channels_list_filepath)"
         << endl;
    return 1;
  }

  if (retrieve_lumi_FCal_) {
    int err_RC = ReadCalibrations(calibrations_filepath_);
    if (err_RC) {
      cerr << "ERROR: in ReadCalibrations(string calibrations_filepath)"
           << endl;
      return 1;
    }
  }

  return 0;
}

int Analysis::ReadCalibrations(string calibrations_filepath) {
// Reads in calibration values for each of the channels being used (those
//   read in with ReadChannelsList)

  for (auto &channel: channels_list_) {
    std::ifstream calibrations_file(calibrations_filepath);
    if (!calibrations_file) {
      cerr << "ERROR: Could not locate calibrations file \'"
           << calibrations_filepath << "\'" << endl;
      return 1;
    }

    string channel_name;
    float slope;
    float intercept;
    bool found_channel = false;
    while (calibrations_file >> channel_name >> slope >> intercept) {
      if (channel.first == channel_name) {
        channel.second.slope = slope;
        channel.second.intercept = intercept;
        found_channel = true;
        break;
      }
    }

    calibrations_file.close();

    if (!found_channel) {
      cerr << "ERROR: Could not locate calibration for channel "
           << channel.first << endl;
      return 2;
    }
  } // Used channels loop

  return 0;
}

int Analysis::ReadChannelsList(string channels_list_filepath) {
// Reads the list of channels to be used in the analysis, and stores them
//   in channels_list_

  std::ifstream channels_list_file(channels_list_filepath);
  if (!channels_list_file) {
    cerr << "ERROR: Could not locate channels list file \'"
         << channels_list_filepath << "\'" << endl;
    return 1;
  }

  string channel_name;
  while ( getline(channels_list_file, channel_name) ) {
    if (channel_name[0] == '#') {
      if (verbose_) cout << "Skipping channel " << channel_name << endl;
      continue;
    }
    Analysis::ChannelCalibration channel_calibration;
    channels_list_.insert(std::make_pair(channel_name, channel_calibration));
  }
  channels_list_file.close();

  return 0;
}

int Analysis::ReadPedestals(string pedestals_dir, string run_name) {
// Reads in pedestal values for each of the channels being used (those
//   read in with ReadChannelsList)

  for (auto &channel: channels_list_) {
    auto pedestals_filepath = pedestals_dir + run_name + ".dat";
    std::ifstream pedestals_file(pedestals_filepath);
    if (!pedestals_file) {
      cerr << "ERROR: Could not locate pedestals file \'"
           << pedestals_filepath << "\'" << endl;
      return 1;
    }

    string channel_name;
    float pedestal;
    float something; // I'm not sure what this value is
    bool found_channel = false;
    while (pedestals_file >> channel_name >> pedestal >> something) {
      if (channel.first == channel_name) {
        channel.second.pedestal = pedestal;
        found_channel = true;
        break;
      }
    }

    pedestals_file.close();

    if (!found_channel) {
      cerr << "ERROR: Could not locate pedestal for run " << run_name
           << ", channel " << channel.first << endl;
      return 2;
    }
  } // Used channels loop

  return 0;
}

void Analysis::ReadParams(string params_filepath) {
// Reads in parameters from json file and assigns their values to member
//   variables

  JSONReader parameter_file(params_filepath);

  verbose_ = parameter_file.get<bool>("verbose");

  // Used to calculate FCal luminosity
  f_rev_ = parameter_file.get<double>("lumi_calculation.f_rev");
  x_sec_ = parameter_file.get<double>("lumi_calculation.x_sec");

  // Value of 0 for reference run means set no corrections
  ref_run_number_ = parameter_file.get<int>("ref_run_number");
  string ref_run_str = std::to_string(ref_run_number_);
  corr_A_ = parameter_file.get<double>("corrections."+ref_run_str+".A");
  corr_C_ = parameter_file.get<double>("corrections."+ref_run_str+".C");
  corr_Avg_ = parameter_file.get<double>("corrections."+ref_run_str+".Avg");

  plot_types_ = parameter_file.get_vector<string>("plot_types");

  // Relevant directories
  calibrations_filepath_ = parameter_file.get<string>("filepaths.calibrations");
  channels_list_filepath_ =
    parameter_file.get<string>("filepaths.channels_list");
  pedestals_dir_ = parameter_file.get<string>("filepaths.pedestals");
  trees_dir_ = parameter_file.get<string>("filepaths.trees");
  run_list_dir_ = parameter_file.get<string>("filepaths.run_list");
  output_dir_ = parameter_file.get<string>("filepaths.output");
}

int Analysis::RunAnalysis() {
// Control flow for the analysis of samples and creation of plots

  ifstream run_names_file(run_list_dir_);
  if (!run_names_file) {
    cerr << "ERROR: Could not locate file: " << run_list_dir_ << endl;
    return 1;
  }

  // Iterate over samples
  string run_name;
  while ( getline(run_names_file, run_name) ) {
    // Skip samples commented out with #
    if (run_name[0] == '#') {
      if (verbose_) cout << "Skipping run " << run_name << endl;
      continue;
    }

    if (retrieve_currents_) {
      if (verbose_) cout << "Retrieving channel pedestals" << endl;
      int err_pedestal = ReadPedestals(pedestals_dir_, run_name);
      if (err_pedestal) {
        cerr << "ERROR: In ReadPedestals(string pedestals_dir, string run_name)"
             << endl;
        return 1;
      }
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
