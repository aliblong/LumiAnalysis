#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TTree.h"

#include "boost/expected/expected.hpp"

#include "analysis.h"
#include "branch_array_buffer_sizes.h"
#include "cutoffs.h"
#include "json_reader.h"
#include "plotter.h"
#include "lumi_current_plot_options.h"
#include "mu_stab_plot_options.h"
#include "single_run_data.h"

using std::string;
using std::vector;

using std::cout;
using std::cerr;
using std::endl;

using boost::expected;
using boost::make_expected;
using boost::make_unexpected;

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

  int nLB = 0;
  tree->SetBranchAddress("NLB", &nLB);
  tree->GetEntry(0);

  return nLB;
}

void DumpVector(const vector<Float_t> &vec) {
  // Prints newline-separated vector elements to cout.
  for (auto element: vec) cout << element << endl;
}

bool IsZeroes(const vector<Float_t> &vec) {
  // Checks whether or not the vector has all elements less than epsilon

  Float_t epsilon = 0.0000000001;
  for (auto element: vec) {
    if (element > epsilon) return false;
  }
  return true;
}

} // Anonymous namespace

Analysis::Analysis(string params_filepath)
  : verbose_(false),
    do_Benedetto_(false),

    f_rev_(0.0),
    x_sec_(0.0),
    ref_run_number_(0),
    corr_A_(0.0),
    corr_C_(0.0),
    corr_Avg_(0.0),

    params_filepath_(""),
    calibrations_filepath_(""),
    channels_list_filepath_(""),
    pedestals_dir_(""),
    trees_dir_(""),
    run_list_dir_(""),
    base_output_dir_(""),
    plots_output_dir_(""),
    benedetto_output_dir_(""),

    retrieve_timestamps_(false),
    retrieve_currents_(false),
    retrieve_lumi_BCM_(false),
    retrieve_lumi_FCal_(false) {
// Initializes various quantities required to run the analysis

  params_filepath_ = params_filepath;
  int err = PrepareAnalysis(params_filepath);
  if (err) {
    cerr << "ERROR: in PrepareAnalysis(string params_filepath)" << endl;
  }
}

int Analysis::AnalyseTree(SingleRunData &this_run) {
// Runs the analysis for a single sample.
// Saves relevant results to member variables of the analysis object.

  string run_name = this_run.run_name_;

  if (verbose_) {
    cout << "Analysing sample " << run_name << endl;
  }

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
  this_run.nLB_ = GetNumLumiBlocks(this_tree);

  // Allocate buffers big enough to hold the data (can't use dynamic memory
  //   because of how TTrees work)
  Float_t currents[gMaxNumChannels][gMaxNumLB];
  Float_t lumi_BCM[gMaxNumLB];
  int quality[gMaxNumLB];
  //std::unique_ptr< vector<string> > beam_mode(nullptr);
  vector<string> *beam_mode = nullptr;

  // Set variables to retrieve based on which plot types are being produced
  if (retrieve_currents_) {
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

  Int_t timestamp; // They saved the timestamp as a signed int on file...
  this_tree->SetBranchAddress("StartOfRun", &timestamp);
  this_tree->SetBranchAddress("quality", &quality);
  this_tree->SetBranchAddress("mode", &beam_mode);
  this_tree->SetBranchAddress("ncoll", &this_run.nCollisions_);

  // Populate those variables which have been set to branches
  this_tree->GetEntry(0);

  this_run.timestamp_ = timestamp;
  if (run_name == "206955") {
    cout << "Manually setting nColl for run 206955" << endl;
    this_run.nCollisions_ = 1368;
  } else if (run_name == "208642") {
    cout << "Manually setting nColl for run 208642" << endl;
    this_run.nCollisions_ = 465;
  } else if (run_name == "211620") {
    cout << "Manually setting nColl for run 211620" << endl;
    this_run.nCollisions_ = 801;
  }

  // Save results to the SingleRunData instance
  for (int iLB = 0; iLB < this_run.nLB_; ++iLB) {
    if (quality[iLB] && beam_mode->at(iLB) == "STABLE BEAMS") {
      if (!this_run.LB_stability_offset_has_been_set_) {
        this_run.LB_stability_offset_ = iLB + 1; //LB numbers are 1-indexed, I think
        this_run.LB_stability_offset_has_been_set_ = true;
        cout << "Timestamp: " << timestamp << endl;
        cout << "First stable LB: " << iLB << endl;
        cout << "Mu for first channel at first stable LB: " << currents[0][iLB]*x_sec_/this_run.nCollisions_/f_rev_ << endl;
      }
      if (retrieve_currents_) {
        unsigned iChannel = 0;
        for (const auto &channel: channels_list_) {
          string this_channel_name = channel.first;

          vector<Float_t> current_vec;
          this_run.currents_.insert(std::make_pair(this_channel_name,
                                                   current_vec));

          Float_t this_channel_pedestal = this_run.pedestals_.at(this_channel_name);
          auto &this_channel_currents = this_run.currents_.at(this_channel_name);

          auto this_current = currents[iChannel][iLB] -
                              this_channel_pedestal;
          this_channel_currents.push_back(this_current);
          ++iChannel;
        }
      }
      if (retrieve_lumi_BCM_) {
        this_run.lumi_BCM_.push_back( lumi_BCM[iLB] );
      }
    }
  }

  //for (const auto &channel: channels_list_) {
  //  auto this_channel_name = channel.first;
  //  auto &this_channel_currents = this_run.currents_.at(this_channel_name);
  //  WriteVectorToFile<Float_t>(this_channel_currents, "output/current_test/"+run_name+"_"+this_channel_name+".root");
  //}
  return 0;
}

int Analysis::CalcFCalLumi(SingleRunData &this_run) {
// Calculates FCal luminosity from currents data for a run
  if (verbose_) cout << "Calculating FCal luminosity" << endl;

  // Checks that these values were not previously calculated
  if (this_run.lumi_FCal_A_.size() > 0 || this_run.lumi_FCal_C_.size() > 0) {
    cerr << "ERROR: FCal lumi vector(s) already filled" << endl;
    return 1;
  }

  // Checks that currents data exists
  if (this_run.currents_.size() == 0) {
    cerr << "ERROR: currents map empty" << endl;
    return 2;
  }
  unsigned nLB_stable = this_run.currents_.begin()->second.size();
  if (nLB_stable == 0) {
    cerr << "ERROR: no stable lumi blocks" << endl;
    return 3;
  }

  // Averages lumi measurement from all channels for each side and for each
  //   lumi block
  for (unsigned iLB = 0; iLB < nLB_stable; ++iLB) {
    Float_t lumi_A_temp = 0.0;
    Float_t lumi_C_temp = 0.0;
    unsigned num_channels_A = 0;
    unsigned num_channels_C = 0;

    // Calculates lumi for a channel using its current and calibration and
    //   adds to running total
    for (const auto &this_module_currents: this_run.currents_) {
      auto current = this_module_currents.second.at(iLB);

      // Skip channel if current is too low
      if (current < gFCalCurrentCutoff) continue;

      string module_name = this_module_currents.first;
      auto intercept = channels_list_.at(module_name).intercept;
      auto slope = channels_list_.at(module_name).slope;

      if (module_name.at(1) == '1') { // module from C-side
        lumi_C_temp += current*slope + intercept;
        ++num_channels_C;
      } else {
        lumi_A_temp += current*slope + intercept;
        ++num_channels_A;
      }
    }

    // Computes the average lumi for this lumi block and adds to lumi vector.
    // Adds 0 if no channels had high enough current
    if (num_channels_A == 0) {
      this_run.lumi_FCal_A_.push_back(0.0);
    } else {
      this_run.lumi_FCal_A_.push_back(corr_A_*lumi_A_temp / num_channels_A);
    }
    if (num_channels_C == 0) {
      this_run.lumi_FCal_C_.push_back(0.0);
    } else {
      this_run.lumi_FCal_C_.push_back(corr_C_*lumi_C_temp / num_channels_C);
    }
  }

  return 0;
}

int Analysis::CalcFCalMu(SingleRunData &this_run) {
// Calculates <mu> from FCal luminosity for a run

  if (verbose_) cout << "Calculating FCal <mu>" << endl;

  // Checks that FCal lumi data exists
  if (this_run.lumi_FCal_A_.size() == 0 && this_run.lumi_FCal_C_.size() == 0) {
    cerr << "ERROR: FCal lumi vector(s) have not been filled" << endl;
    return 1;
  }
  // Checks that mu data has not already been calculated
  if (this_run.mu_FCal_A_.size() > 0 || this_run.mu_FCal_C_.size() > 0) {
    cerr << "ERROR: FCal mu vector(s) already filled" << endl;
    return 2;
  }
  // Checks that value for number of collisions is nonzero
  if (this_run.nCollisions_ == 0) {
    cerr << "ERROR: nCollisions = 0" << endl;
    return 3;
  }

  Float_t conversion_factor = x_sec_ / (this_run.nCollisions_ * f_rev_);
  // Cross-section conversion factor from 7 -> 8 TeV
  if (this_run.run_name_.at(0) == '2') conversion_factor /= 1.05;

  // Calculates <mu> for each lumi block
  for (const auto &lumi: this_run.lumi_FCal_A_) {
    this_run.mu_FCal_A_.push_back(lumi*conversion_factor);
  }
  for (const auto &lumi: this_run.lumi_FCal_C_) {
    this_run.mu_FCal_C_.push_back(lumi*conversion_factor);
  }
  return 0;
}

int Analysis::CreateAllRunPlots(const std::map<string, SingleRunData> &runs_data) {
  for (const auto &plot_type: plot_types_) {
    if (plot_type == "mu_stability") {
      if (verbose_) cout << "Making mu stability plot" << endl;
      MuStabPlotOptions plot_options(params_filepath_);
      Plotter::PlotMuStability(runs_data, plot_options, plots_output_dir_);
    }
  }

  return 0;
}

int Analysis::CreateBenedettoOutput(const SingleRunData &this_run) const {
  int err = this_run.CreateBenedettoOutput(benedetto_output_dir_);
  return err;
}

int Analysis::CreateLumiCurrentPlots(const SingleRunData &this_run) {
  if (verbose_) cout << "    " << "Making lumi vs. current plots" << endl;

  string run_name = this_run.run_name_;
  LumiCurrentPlotOptions plot_options(params_filepath_);
  std::map<string, FitResults> fit_results;

  if (plot_options.do_individual()) {
    for (const auto &channel: channels_list_) {
      auto this_channel_name = channel.first;
      auto this_channel_current = this_run.currents_.at(this_channel_name);

      if (IsZeroes(this_channel_current) && verbose_) {
        if (verbose_) cout << "Skipping channel with zero current: "
                           << this_channel_name << endl;
        continue;
      }

      FitResults this_channel_fit_results;
      if (this_run.pedestals_.at(this_channel_name) > 20) {
        this_channel_fit_results.is_short = true;
      } else {
        this_channel_fit_results.is_short = false;
      }

      auto err_plot = Plotter::PlotLumiCurrent(this_run.lumi_BCM_,
                                              this_channel_current,
                                              run_name,
                                              this_channel_name,
                                              plot_options,
                                              plots_output_dir_,
                                              this_channel_fit_results);
      //TODO
      //if (err_plot)

      if (plot_options.do_fit()) {
        fit_results.insert(std::make_pair(this_channel_name,
                                          this_channel_fit_results));
      }
    } //channels loop
  }

  if (plot_options.do_sum()) {
    vector<Float_t> channel_currents_sum_A;
    vector<Float_t> channel_currents_sum_C;
    for (const auto &channel: channels_list_) {
      auto this_channel_name = channel.first;
      auto this_channel_current = this_run.currents_.at(this_channel_name);

      vector<Float_t> *this_side_channel_currents_sum;
      if (this_channel_name.at(1) == '1') {
        this_side_channel_currents_sum = &channel_currents_sum_A;
      } else if (this_channel_name.at(1) == '8') {
        this_side_channel_currents_sum = &channel_currents_sum_C;
      } else {
        cerr << "Invalid FCal module in Analysis::CreateSingleRunPlots" << endl;
        return 1;
      }
      if (this_side_channel_currents_sum->size() == 0) {
        *this_side_channel_currents_sum = this_channel_current;
      } else {
        std::transform(this_side_channel_currents_sum->begin(),
                       this_side_channel_currents_sum->end(),
                       this_channel_current.begin(),
                       this_side_channel_currents_sum->begin(),
                       std::plus<Float_t>());
      }
    }

    FitResults channels_sum_fit_results;
    /*
    auto err_plot = Plotter::PlotLumiTotalCurrent(
                              this_run.lumi_BCM_,
                              channel_currents_sum_A,
                              channel_currents_sum_C,
                              run_name,
                              plot_options,
                              plots_output_dir_);
                                            */
    auto err_plot = Plotter::PlotLumiCurrent(this_run.lumi_BCM_,
                                            channel_currents_sum_A,
                                            run_name,
                                            "Sum_A",
                                            plot_options,
                                            plots_output_dir_,
                                            channels_sum_fit_results);
    err_plot = Plotter::PlotLumiCurrent(this_run.lumi_BCM_,
                                            channel_currents_sum_C,
                                            run_name,
                                            "Sum_C",
                                            plot_options,
                                            plots_output_dir_,
                                            channels_sum_fit_results);
  }

  if (plot_options.do_fit()) {
    auto err_fit_res = Plotter::WriteFitResultsToTree(fit_results,
                                                     run_name,
                                                     fit_results_output_dir_);

    auto err_calib = Plotter::WriteCalibrationToText(fit_results,
                                                    run_name,
                                                    calibrations_output_dir_);

    auto result_geo = Plotter::GeometricAnalysisOfFitResults(fit_results,
                                                             run_name,
                                                             geometric_analysis_output_dir_);

    if (!result_geo.valid()) cout << result_geo.error()->what() << endl;
  }

  return 0;
}

int Analysis::CreateSingleRunPlots(const SingleRunData &this_run) {
// Create those plots which use data from only a single sample

  string run_name = this_run.run_name_;

  if (plot_types_.size() == 0) {
    if (verbose_) cout << "No plot types selected" << endl;
    return 0;
  }
  if (verbose_) cout << "Creating plots for sample " << run_name << endl;

  for (const auto &plot_type: plot_types_) {
    if (plot_type == "lumi_current") {
      std::map<string, FitResults> fit_results;
      int err_lumi_curr = CreateLumiCurrentPlots(this_run);
    }
  }

  return 0;
}

int Analysis::PrepareAnalysis(string params_filepath) {
// Reads in the analysis parameters (json) and channel information (text). 
// Sets flags for which data to retrieve based on which plots to produce.

  ReadParams(params_filepath);
  for (const auto &plot_type: plot_types_) {
    if (plot_type == "lumi_current" ) {
      retrieve_currents_ = true;
      retrieve_lumi_BCM_ = true;
    }
    else if (plot_type == "mu_stability") {
      retrieve_timestamps_ = true;
      retrieve_lumi_BCM_ = true;
      retrieve_lumi_FCal_ = true;
    }
  }

  if (do_Benedetto_) {
    retrieve_lumi_FCal_ = true;
  }

  // FCal currents are required to calculate FCal lumi
  if (retrieve_lumi_FCal_) retrieve_currents_ = true;

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
    Float_t slope;
    Float_t intercept;
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

void Analysis::ReadParams(string params_filepath) {
// Reads in parameters from json file and assigns their values to member
//   variables

  JSONReader parameter_file(params_filepath);

  verbose_ = parameter_file.get<bool>("verbose");

  // Text output file for Benedetto
  do_Benedetto_ = parameter_file.get<bool>("do_Benedetto");

  // Used to calculate FCal luminosity
  //   Bunch crossing frequency
  f_rev_ = parameter_file.get<double>("lumi_calculation.f_rev");
  //   Cross-section for pp interaction @ 7 TeV
  x_sec_ = parameter_file.get<double>("lumi_calculation.x_sec");

  // Value of 0 for reference run means set no corrections
  ref_run_number_ = parameter_file.get<int>("ref_run_number");
  string ref_run_str = std::to_string(ref_run_number_);
  corr_A_ = parameter_file.get<double>("corrections."+ref_run_str+".A");
  corr_C_ = parameter_file.get<double>("corrections."+ref_run_str+".C");
  corr_Avg_ = parameter_file.get<double>("corrections."+ref_run_str+".Avg");

  auto plot_types_map = parameter_file.get_map<bool>("plot_types");
  for (const auto &plot_type: plot_types_map) {
    if (plot_type.second == true) plot_types_.push_back(plot_type.first);
  }

  // Relevant directories
  calibrations_filepath_ = parameter_file.get<string>("filepaths.calibrations");
  channels_list_filepath_ =
    parameter_file.get<string>("filepaths.channels_list");
  pedestals_dir_ = parameter_file.get<string>("filepaths.pedestals");
  trees_dir_ = parameter_file.get<string>("filepaths.trees");
  run_list_dir_ = parameter_file.get<string>("filepaths.run_list");
  base_output_dir_ = parameter_file.get<string>("filepaths.output.base");
  plots_output_dir_ = base_output_dir_ +
                     parameter_file.get<string>("filepaths.output.plots");
  benedetto_output_dir_ = base_output_dir_ +
                          parameter_file.get<string>(
                            "filepaths.output.benedetto");

  fit_results_output_dir_ = base_output_dir_ +
                            parameter_file.get<string>(
                              "filepaths.output.fit_results.base");
  calibrations_output_dir_ = fit_results_output_dir_ +
                             parameter_file.get<string>(
                               "filepaths.output.fit_results.calibrations");
  geometric_analysis_output_dir_ = fit_results_output_dir_ +
                                   parameter_file.get<string>(
                                     "filepaths.output.fit_results.geometric");
}

int Analysis::RunAnalysis() {
// Control flow for the analysis of samples and creation of plots

  std::ifstream run_names_file(run_list_dir_);
  if (!run_names_file) {
    cerr << "ERROR: Could not locate file: " << run_list_dir_ << endl;
    return 1;
  }

  string run_name;
  vector<string> run_names_vec;
  while ( getline(run_names_file, run_name) ) {
    // Skip samples commented out with #
    if (run_name[0] == '#') {
      if (verbose_) cout << "Skipping run " << run_name << endl;
      continue;
    }
    run_names_vec.push_back(run_name);
  }

  vector<string> channel_names;
  for (const auto &channel: channels_list_) {
    channel_names.push_back(channel.first);
  }


  // Iterate over samples
  std::map<string, SingleRunData> runs_data;

  //#pragma omp parallel for
  for (auto it = run_names_vec.begin(); it < run_names_vec.end(); ++it) {
    SingleRunData this_run(*it);

    // Read in pedestals
    if (retrieve_currents_) {
      if (verbose_) cout << "Retrieving channel pedestals" << endl;

      int err_pedestal = this_run.ReadPedestals(pedestals_dir_, channel_names);
      if (err_pedestal) {
        cerr << "ERROR: In ReadPedestals(string pedestals_dir, string run_name)"
             << endl;
        continue;
        //return 1;
      }
    }

    auto err_analyse = AnalyseTree(this_run);
    if (err_analyse) {
      cerr << "ERROR: In AnalyseTree(string run_name)" << endl;
      continue;
      //return 1;
    }

    auto err_plots = CreateSingleRunPlots(this_run);
    if (err_plots) {
      cerr << "ERROR: In CreateSingleRunPlots(string run_name)" << endl;
      continue;
      //return 1;
    }

    if (retrieve_lumi_FCal_) {
      auto err_lumi_FCal = CalcFCalLumi(this_run);
      auto err_mu_FCal = CalcFCalMu(this_run);
    }
    if (do_Benedetto_) {
      auto err_Benedetto = CreateBenedettoOutput(this_run);
    }

    //#pragma omp critical
    runs_data.insert(std::make_pair(*it, std::move(this_run)));
  }

  auto err_plots_all = CreateAllRunPlots(runs_data);

  return 0;
}

/*
int WriteCurrentsToFile(vector<Float_t> currents, string run_name) {
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
