#ifndef INCLUDE_FCAL_REGION_Z_SIDE_H_
#define INCLUDE_FCAL_REGION_Z_SIDE_H_

#include <set>

namespace FCalRegion {

enum class ZSide {A, C, Both};

const std::set<ZSide> Z_SIDES = {ZSide::A,
                                 ZSide::C,
                                 ZSide::Both};
}

#endif

