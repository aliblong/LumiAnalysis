#ifndef LUMIANALYSIS_INCLUDE_SINGLE_RUN_DATA_H_
#define LUMIANALYSIS_INCLUDE_SINGLE_RUN_DATA_H_

#include <string>
#include <vector>

#include "boost/container/flat_map.hpp"

#include "Rtypes.h"

#include "error.h"
#include "void.h"

class Analysis;

class SingleRunData {
 public:
  //friend class Analysis;

  SingleRunData(std::string run_name, const Analysis* analysis);
  ~SingleRunData() {}

  Error::Expected<Void> Init();

  Error::Expected<Void> CreateLumiCurrentPlots() const;
  Error::Expected<Void> CreateSingleRunPlots() const;
  Error::Expected<Void> CalcFCalLumi();
  Error::Expected<Void> CalcFCalMu();
  Error::Expected<Void> CreateBenedettoOutput() const;

  auto timestamp() const { return timestamp_; }
  auto run_name() const { return run_name_; }
  auto nCollisions() const { return nCollisions_; }
  auto nLB() const { return nLB_; }
  auto LB_stability_offset() const { return LB_stability_offset_; }

  const auto& pedestals() const { return pedestals_; }
  const auto& currents() const { return currents_; }
  const auto& lumi_ofl() const { return lumi_ofl_; }
  const auto& lumi_FCal_A() const { return lumi_FCal_A_; }
  const auto& lumi_FCal_C() const { return lumi_FCal_C_; }
  const auto& mu_FCal_A() const { return mu_FCal_A_; }
  const auto& mu_FCal_C() const { return mu_FCal_C_; }

 private:
  Error::Expected<Void> ReadPedestals();
  Error::Expected<Void> ReadTree();
  void InitCurrentsMap();
  void HardcodenLBIfMissingFromTree();

  std::string run_name_;
  // This reference allows access to analysis-wide parameters such as output
  //   directories.
  const Analysis* analysis_;

  Int_t nLB_= 0; // Number of lumi blocks
  Int_t nCollisions_ = 0;
  Int_t timestamp_ = 0;

  Int_t LB_stability_offset_ = 0;

  boost::container::flat_map<std::string, Float_t> pedestals_;
  boost::container::flat_map< std::string, std::vector<Float_t> > currents_;
  std::vector<Float_t> lumi_ofl_;
  std::vector<Float_t> lumi_FCal_A_;
  std::vector<Float_t> lumi_FCal_C_;
  std::vector<Float_t> mu_FCal_A_;
  std::vector<Float_t> mu_FCal_C_;
};

#endif
