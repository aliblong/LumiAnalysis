#ifndef LUMIANALYSIS_INCLUDE_FCAL_REGION_DATA_H_
#define LUMIANALYSIS_INCLUDE_FCAL_REGION_DATA_H_

#include <string>
#include <memory>

#include "fcal_region_z_side.h"
#include "TProfile.h"

class FCalRegionData {
 public:
  FCalRegionData(FCalRegion::ZSide region,
                 std::string&& name,
                 std::string&& title,
                 unsigned color,
                 Float_t size,
                 unsigned style);
  ~FCalRegionData();
  FCalRegionData(const FCalRegionData &rhs);
  FCalRegionData(FCalRegionData &&rhs);

  FCalRegion::ZSide region_name_;

  std::string plot_name_;
  std::string plot_title_;

  std::unique_ptr<TProfile> plot_;

  unsigned marker_color_ = 0;
  Float_t marker_size_ = 1.0;
  unsigned marker_style_ = 0;
};

#endif
