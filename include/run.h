#ifndef LUMIANALYSIS_INCLUDE_SINGLE_RUN_DATA_H_
#define LUMIANALYSIS_INCLUDE_SINGLE_RUN_DATA_H_

#include <string>
#include <vector>

#include "boost/container/flat_map.hpp"

#include "Rtypes.h"

#include "analysis.h"
#include "error.h"
#include "void.h"

class Run {
 public:
  Run(std::string run_name, const Analysis* analysis);
  ~Run() {}

  Error::Expected<Void> Init();

  // Create all plots which depend on data from a single run only
  Error::Expected<Void> CreateSingleRunPlots() const;
  // Current vs. some luminosity algorithm from which we can derive a calibration for 
  //   current->luminosity
  Error::Expected<Void> CreateLumiCurrentPlots() const;
  // FCal luminosity and mu values are calculated by applying a derived calibration to the current
  //   values
  Error::Expected<Void> CalcFCalLumi();
  Error::Expected<Void> CalcFCalMu();
  // Benedetto output, aka BG files, list for each lumi block on a new line the lumi block, A-side
  //   mu, and C-side mu
  Error::Expected<Void> CreateBenedettoOutput() const;

  // Timestamp for the start of the run
  auto timestamp() const { return timestamp_; }
  // I.e. run number
  auto run_name() const { return run_name_; }
  // Number of colliding bunches in the LHC
  auto nBunches() const { return nBunches_; }
  // Number of luminosity blocks - arbitrary divisions of time typically equal to 1 min
  auto nLB() const { return nLB_; }
  // First ready-for-physics (RFP) LB number
  auto LB_stability_offset() const { return LB_stability_offset_; }
  // Beamspot z position averaged over the entire run
  float avg_beamspot_z() const { return avg_beamspot_z_; }
  // Default val for beamspot in LB where it is absent
  static auto BeamspotPlaceholderVal() { return -9999.0; }

  // Vectors containing a value for each LB
  const auto& currents() const { return currents_; }
  const auto& lumi_ofl() const { return lumi_ofl_; }
  const auto& lumi_FCal_A() const { return lumi_FCal_A_; }
  const auto& lumi_FCal_C() const { return lumi_FCal_C_; }
  const auto& mu_ofl() const { return mu_ofl_; }
  const auto& mu_FCal_A() const { return mu_FCal_A_; }
  const auto& mu_FCal_C() const { return mu_FCal_C_; }
  const auto& beamspot_z() const { return beamspot_z_; }

  // Maps over channels
  const auto& pedestals() const { return pedestals_; }
  const auto& channel_calibrations() const { return channel_calibrations_; }

 private:
  Error::Expected<Void> ReadPedestals();
  Error::Expected<Void> ReadTree();
  // Determine which LB range to use based on user input first then RFP flag
  Error::Expected<std::array<Int_t,2>> GetLBBounds() const;
  // Read in nBunches from an external file; this value is normally retrieved and stored during
  //   tree creation, but a bug results in the retrieved value being 0 in ~10% of runs
  void GetExternalNBunches();

  std::string run_name_;
  // This allows access to analysis-wide parameters e.g. output directories.
  const Analysis* analysis_;

  Int_t nLB_= 0;
  Int_t nBunches_ = 0;
  Int_t timestamp_ = 0;

  Int_t LB_stability_offset_ = 0;

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
