#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "single_run_data.h"

using std::string;
using std::vector;

using std::cout;
using std::cerr;
using std::endl;

int SingleRunData::ReadPedestals(string pedestals_dir,
                                 const vector<string> &channel_names) {
// Reads in pedestal values for each of the channels being used (those
//   read in with ReadChannelsList)

  for (const auto &this_channel_name: channel_names) {
    auto pedestals_filepath = pedestals_dir + run_name_ + ".dat";
    std::ifstream pedestals_file(pedestals_filepath);
    if (!pedestals_file) {
      cerr << "ERROR: Could not locate pedestals file \'"
           << pedestals_filepath << "\'" << endl;
      return 1;
    }

    string channel_name;
    Float_t pedestal;
    Float_t something; // I'm not sure what this value is
    bool found_channel = false;
    while (pedestals_file >> channel_name >> pedestal >> something) {
      if (this_channel_name == channel_name) {
        pedestals_.insert(std::make_pair(channel_name, pedestal));
        found_channel = true;
        break;
      }
    }

    pedestals_file.close();

    if (!found_channel) {
      cerr << "ERROR: Could not locate pedestal for run " << run_name_
           << ", channel " << this_channel_name << endl;
      return 2;
    }
  } // Used channels loop

  return 0;
}

int SingleRunData::CreateBenedettoOutput(string output_dir) const {
  int err_system = system( ("mkdir -p "+output_dir).c_str() );
  auto out_filepath = output_dir+run_name_+".dat";
  std::ofstream out_file(out_filepath);
  if (!out_file.is_open()) {
    cerr << "ERROR: Could not open file \'" << out_filepath << "\'" << endl;
    return 1;
  }

  auto num_points = mu_FCal_A_.size();
  for (unsigned iLB = 0; iLB < num_points; ++iLB) {
    out_file << (iLB + LB_stability_offset_) << ' ' << mu_FCal_A_.at(iLB) << ' '
             << mu_FCal_C_.at(iLB) << std::endl;
  }
  out_file.close();
  return 0;
}
