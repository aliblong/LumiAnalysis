#ifndef INCLUDE_FCAL_REGION_Z_SIDE_H_
#define INCLUDE_FCAL_REGION_Z_SIDE_H_

#include "enum.h"

namespace FCalRegion {

class ZSide: public Enum<ZSide> {
 public:
  using Enum::Enum;
  //explicit ZSide(int val) : Enum(val) {}

  static const ZSide A;
  static const ZSide C;
  static const ZSide Both;
};

}

#endif

