#ifndef INCLUDE_FCAL_REGION_H_
#define INCLUDE_FCAL_REGION_H_

#include <set>
#include <string>

#include "Rtypes.h"

#include "boost/expected/expected.hpp"

#include "error.h"
#include "fcal_region_z_side.h"

namespace FCalRegion {

enum class Axis {X, Y};

enum class Sign {Pos, Neg};

//bool operator< (const Axis &a1, const Axis &a2) {
//  if (a2 == Axis::X) return false;
//  else {
//    if (a1 == Axis::Y) return false;
//    else return true;
//  }
//}

//bool operator< (const Sign &s1, const Sign &s2) {
//  if (s2 == Sign::Pos) return false;
//  else {
//    if (s1 == Sign::Neg) return false;
//    else return true;
//  }
//}

const std::set<Axis> AXES = {Axis::X, Axis::Y};
const std::set<Sign> SIGNS = {Sign::Pos, Sign::Neg};

typedef std::tuple<ZSide, Axis, Sign> ModuleHalf;

typedef std::set<ModuleHalf> ModuleHalfSet;
ModuleHalfSet CreateModuleHalfSet();

std::string ToString(ZSide region);
std::string ToString(Axis axis);
std::string ToString(Sign sign);

Error::Expected<std::string> PhiSliceFromChannel(std::string channel_name);
Sign SignFromAxisAndPhiSlice(Axis axis, std::string phi_slice_name);

}

#endif

//void DumpGeoRegion(const ModuleHalf &region);
