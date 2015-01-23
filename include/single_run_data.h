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
  friend class Analysis;

  SingleRunData(std::string run_name)
      : run_name_(run_name),
        nLB_(0),
        nCollisions_(0),
        LB_stability_offset_(0),
        LB_stability_offset_has_been_set_(false) {};
  ~SingleRunData(){};

  auto timestamp() const { return timestamp_; }
  auto run_name() const { return run_name_; }
  auto nCollisions() const { return nCollisions_; }
  auto nLB() const { return nLB_; }
  auto LB_stability_offset() const { return LB_stability_offset_; }
  auto LB_stability_offset_has_been_set() const { return LB_stability_offset_has_been_set_; }

  const auto& pedestals() const { return pedestals_; }
  const auto& currents() const { return currents_; }
  const auto& lumi_BCM() const { return lumi_BCM_; }
  const auto& lumi_FCal_A() const { return lumi_FCal_A_; }
  const auto& lumi_FCal_C() const { return lumi_FCal_C_; }
  const auto& mu_FCal_A() const { return mu_FCal_A_; }
  const auto& mu_FCal_C() const { return mu_FCal_C_; }

 private:
  Error::Expected<Void> ReadPedestals(std::string pedestals_dir,
                    const std::vector<std::string> &channel_names);
  Error::Expected<Void> CreateBenedettoOutput(std::string output_dir) const;

  std::string run_name_;
  Int_t nLB_; // Number of lumi blocks
  Int_t nCollisions_;
  Long64_t timestamp_;

  Int_t LB_stability_offset_;
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
