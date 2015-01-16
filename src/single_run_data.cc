#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "boost/expected/expected.hpp"

#include "error.h"
#include "util.h"
#include "single_run_data.h"

using std::string;
using std::vector;

using std::cout;
using std::cerr;
using std::endl;

using std::make_shared;

using boost::make_unexpected;

using Error::Expected;

// Reads in pedestal values for each of the channels being used (those
//   read in with ReadChannelsList)
Expected<Void> SingleRunData::ReadPedestals(string pedestals_dir,
                                 const vector<string> &channel_names) {
  auto this_func_name = "SingleRunData::ReadPedestals";

  for (const auto &this_channel_name: channel_names) {
    auto pedestals_filepath = pedestals_dir + run_name_ + ".dat";
    std::ifstream pedestals_file(pedestals_filepath);
    if (!pedestals_file) {
      return make_unexpected(make_shared<Error::File>(pedestals_filepath, this_func_name));
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
      auto err_msg = "could not locate pedestal for run "+run_name_;
      return make_unexpected(make_shared<Error::Runtime>(err_msg, this_func_name));
    }
  }

  return Void();
}

Expected<Void> SingleRunData::CreateBenedettoOutput(string output_dir) const {
  auto this_func_name = "SingleRunData::CreateBenedettoOutput";

  TRY( Util::mkdir(output_dir) )
  auto out_filepath = output_dir+run_name_+".dat";
  std::ofstream out_file(out_filepath);
  if (!out_file.is_open()) {
    return make_unexpected(make_shared<Error::File>(out_filepath, this_func_name));
  }

  auto num_points = mu_FCal_A_.size();
  for (unsigned iLB = 0; iLB < num_points; ++iLB) {
    out_file << (iLB + LB_stability_offset_) << ' ' << mu_FCal_A_.at(iLB) << ' '
             << mu_FCal_C_.at(iLB) << std::endl;
  }
  out_file.close();
  return Void();
}
