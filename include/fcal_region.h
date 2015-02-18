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

const std::set<Axis> AXES = {Axis::X, Axis::Y};
const std::set<Sign> SIGNS = {Sign::Pos, Sign::Neg};

typedef std::tuple<ZSide, Axis, Sign> ModuleHalf;

typedef std::set<ModuleHalf> ModuleHalfSet;
ModuleHalfSet CreateModuleHalfSet();

std::string ToString(ZSide region);
std::string ToString(Axis axis);
std::string ToString(Sign sign);

Error::Expected<std::string> PhiSliceFromChannel(const std::string& channel_name);
Sign SignFromAxisAndPhiSlice(Axis axis, const std::string& phi_slice_name);

bool IsValidChannel(const std::string& channel_name);

ZSide ToZSide(const std::string& channel_name);

}

#endif

//void DumpGeoRegion(const ModuleHalf &region);
