#include "fcal_region_data.h"
#include <string>
#include "Rtypes.h"

FCalRegionData::FCalRegionData(FCalRegion::ZSide region,
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
      marker_style_(style),
      plot_(nullptr) {}

FCalRegionData::~FCalRegionData() {}

FCalRegionData::FCalRegionData(const FCalRegionData &rhs)
    : region_name_(rhs.region_name_),
      plot_name_(rhs.plot_name_),
      plot_title_(rhs.plot_title_),
      marker_color_(rhs.marker_color_),
      marker_size_(rhs.marker_size_),
      marker_style_(rhs.marker_style_),
      plot_(nullptr) {}

FCalRegionData::FCalRegionData(FCalRegionData &&rhs)
    : region_name_(rhs.region_name_),
      plot_name_(rhs.plot_name_),
      plot_title_(rhs.plot_title_),
      marker_color_(rhs.marker_color_),
      marker_size_(rhs.marker_size_),
      marker_style_(rhs.marker_style_),
      plot_(std::move(rhs.plot_)) {}
