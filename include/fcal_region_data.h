#ifndef LUMIANALYSIS_INCLUDE_FCAL_REGION_DATA_H_
#define LUMIANALYSIS_INCLUDE_FCAL_REGION_DATA_H_

#include <string>

#include "TProfile.h"

enum class FCalRegion {A, C, Avg};

class FCalRegionData {
 public:
  FCalRegionData(FCalRegion region,
                 std::string name,
                 std::string title,
                 unsigned color,
                 Float_t size,
                 unsigned style)
    : region_name_(region),
      plot_name_(name),
      plot_title_(title),
      marker_color_(color),
      marker_size_(size),
      marker_style_(style) {}
  ~FCalRegionData() {}

  FCalRegion region_name_;

  std::string plot_name_;
  std::string plot_title_;

  std::unique_ptr<TProfile> plot_;

  unsigned marker_color_;
  Float_t marker_size_;
  unsigned marker_style_;
};

#endif
