#ifndef LUMIANALYSIS_INCLUDE_SINGLE_RUN_DATA_H_
#define LUMIANALYSIS_INCLUDE_SINGLE_RUN_DATA_H_

#include <string>
#include <vector>

#include "boost/container/flat_map.hpp"

#include "Rtypes.h"

#include "analysis.h"
#include "error.h"
#include "void.h"

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
  auto nBunches() const { return nBunches_; }
  auto nLB() const { return nLB_; }
  auto LB_stability_offset() const { return LB_stability_offset_; }
  float avg_beamspot_z() const { return avg_beamspot_z_; }
  static auto BeamspotPlaceholderVal() { return -9999.0; }
  //Default val for beamspot in LB where it is absent

  const auto& pedestals() const { return pedestals_; }
  const auto& currents() const { return currents_; }
  const auto& lumi_ofl() const { return lumi_ofl_; }
  const auto& lumi_FCal_A() const { return lumi_FCal_A_; }
  const auto& lumi_FCal_C() const { return lumi_FCal_C_; }
  const auto& mu_ofl() const { return mu_ofl_; }
  const auto& mu_FCal_A() const { return mu_FCal_A_; }
  const auto& mu_FCal_C() const { return mu_FCal_C_; }
  const auto& channel_calibrations() const { return channel_calibrations_; }
  const auto& beamspot_z() const { return beamspot_z_; }

 private:
  Error::Expected<Void> ReadPedestals();
  Error::Expected<Void> ReadTree();
  Error::Expected<std::array<Int_t,2>> GetLBBounds() const;
  void InitCurrentsMap();
  void GetExternalNBunches();

  std::string run_name_;
  // This allows access to analysis-wide parameters such as output directories.
  const Analysis* analysis_;

  Int_t nLB_= 0; // Number of lumi blocks
  Int_t nBunches_ = 0; // Number of colliding bunches
  Int_t timestamp_ = 0;

  Int_t LB_stability_offset_ = 0; // First ready-for-physics (RFP) LB number

  boost::container::flat_map<std::string, Float_t> pedestals_;
  boost::container::flat_map< std::string, std::vector<Float_t> > currents_;
  std::vector<Float_t> lumi_ofl_;
  std::vector<Float_t> lumi_FCal_A_;
  std::vector<Float_t> lumi_FCal_C_;
  std::vector<Float_t> mu_ofl_;
  std::vector<Float_t> mu_FCal_A_;
  std::vector<Float_t> mu_FCal_C_;
  std::vector<Int_t> RFP_flag_;
  std::vector<Float_t> beamspot_z_;
  Float_t avg_beamspot_z_ = -999.0;
  boost::container::flat_map<std::string, Analysis::ChannelCalibration> channel_calibrations_;
};

#endif
