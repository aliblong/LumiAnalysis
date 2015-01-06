#ifndef INCLUDE_FCAL_REGION_H_
#define INCLUDE_FCAL_REGION_H_

#include <set>
#include <string>

#include "Rtypes.h"

#include "enum.h"
#include "fcal_region_z_side.h"

namespace FCalRegion {

class Axis: public Enum<Axis> {
 public:
  using Enum::Enum;
  //explicit Axis(int val) : Enum(val) {}

  static const Axis X;
  static const Axis Y;
};

class Sign: public Enum<Sign> {
 public:
  using Enum::Enum;
  //explicit Sign(int val) : Enum(val) {}

  static const Sign Pos;
  static const Sign Neg;
};

typedef std::tuple<ZSide, Axis, Sign> ModuleHalf;

struct ModuleHalfOrderer {
  bool operator() (const ModuleHalf &r1, const ModuleHalf &r2) const {
    auto side1 = std::get<0>(r1);
    auto side2 = std::get<0>(r2);
    auto axis1 = std::get<1>(r1);
    auto axis2 = std::get<1>(r2);
    auto sign1 = std::get<2>(r1);
    auto sign2 = std::get<2>(r2);

    if (side1 == side2) {
      if (sign1 == sign2) {
        return axis1.value() < axis2.value();
      } else {
        return sign1.value() < sign2.value();
      }
    } else {
      return side1.value() < side2.value();
    }
  }
};

typedef std::set<ModuleHalf, ModuleHalfOrderer> ModuleHalfSet;
ModuleHalfSet CreateModuleHalfSet();

std::string ToString(ZSide region);
std::string ToString(Axis axis);
std::string ToString(Sign sign);

std::string PhiSliceFromChannel(std::string channel_name);
Sign SignFromAxisAndPhiSlice(Axis axis, std::string phi_slice_name);

}

#endif

//void DumpGeoRegion(const ModuleHalf &region);
