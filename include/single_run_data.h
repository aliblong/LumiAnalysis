#ifndef LUMIANALYSIS_INCLUDE_SINGLE_RUN_DATA_H_
#define LUMIANALYSIS_INCLUDE_SINGLE_RUN_DATA_H_

#include <string>
#include <vector>
#include <map>

#include "Rtypes.h"

#include "error.h"
#include "void.h"

// Holds data for single runs
class SingleRunData {
 public:
  SingleRunData(std::string run_name)
      : run_name_(run_name),
        nLB_(0),
        nCollisions_(0),
        LB_stability_offset_(0),
        LB_stability_offset_has_been_set_(false) {};
  ~SingleRunData(){};

  ULong_t timestamp() const { return timestamp_; }

 //private:
  friend class Analysis;

  Error::Expected<Void> ReadPedestals(std::string pedestals_dir,
                    const std::vector<std::string> &channel_names);
  Error::Expected<Void> CreateBenedettoOutput(std::string output_dir) const;

  std::string run_name_;
  int nLB_; //number of lumi blocks
  int nCollisions_;
  UInt_t timestamp_;

  int LB_stability_offset_;
  bool LB_stability_offset_has_been_set_;

  std::map<std::string, Float_t> pedestals_;
  std::map< std::string, std::vector<Float_t> > currents_;
  std::vector<Float_t> lumi_BCM_;
  std::vector<Float_t> lumi_FCal_A_;
  std::vector<Float_t> lumi_FCal_C_;
  std::vector<Float_t> mu_FCal_A_;
  std::vector<Float_t> mu_FCal_C_;
};

#endif
