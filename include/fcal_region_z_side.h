#ifndef INCLUDE_FCAL_REGION_Z_SIDE_H_
#define INCLUDE_FCAL_REGION_Z_SIDE_H_

#include <set>

namespace FCalRegion {

enum class ZSide {A, C, Both};

//bool operator< (const ZSide &s1, const ZSide &s2) {
//  if (s2 == ZSide::A) return false;
//  else if (s2 == ZSide::C) {
//    if (s1 == ZSide::A) return true;
//    else return false;
////  }
////  else {
//    if (s1 == ZSide::Both) return false;
//    else return true;
//  }
//}
//
const std::set<ZSide> Z_SIDES = {ZSide::A,
                                 ZSide::C,
                                 ZSide::Both};
}

#endif

