#include <iostream>
#include <map>
#include <sstream>
#include <string>

#include "fcal_region.h"
#include "enum.h"
#include "error.h"

using std::string;

using FCalRegion::ZSide;
using FCalRegion::Axis;
using FCalRegion::Sign;

namespace {

template <typename T>
std::string MsgFromDict(const T &key, std::map<T, std::string> dict) {
  auto msg_it = dict.find(key);
  if (msg_it == dict.end()) {
    return "ERROR: entry for enum not found";
  } else {
    return msg_it->second;
  }
}

string StringFromInt(int i, unsigned width) {
  std::stringstream result;
  result.width(2);
  result.fill('0');
  result << i;
  return result.str();
}

int IntFromChar(char c) {
  int result = c - '0';
  return result;
}

string PhiSliceFromID(char region_ID, char module_ID, char channel_ID) {
  auto module_ID_int = IntFromChar(module_ID);
  if ( (module_ID_int < 0) || (module_ID_int > 9) ) {
    Error::Report("ERROR: in IntFromChar(char c) - module_ID input falls"
                  "outside [0,9]");
    return "ERROR";
  }
  auto channel_ID_int = IntFromChar(channel_ID);
  if ( (channel_ID_int < 0) || (channel_ID_int > 9) ) {
    Error::Report("ERROR: in IntFromChar(char c) - channel_ID input falls"
                  "outside [0,9]");
    return "ERROR";
  }
  unsigned phi_slice_ID_int = (module_ID_int * 2) + (channel_ID_int / 4);
  string region_ID_str(1,region_ID);
  string phi_slice = region_ID_str + "1." + StringFromInt(phi_slice_ID_int, 2);
  return phi_slice;
}

}

INIT_ENUM(Axis);

const Axis Axis::X(0);
const Axis Axis::Y(1);

INIT_ENUM(Sign);

const Sign Sign::Pos(0);
const Sign Sign::Neg(1);

FCalRegion::ModuleHalfSet FCalRegion::CreateModuleHalfSet() {
  FCalRegion::ModuleHalfSet result;
  ENUM_FOR_EACH(side_ptr_it, FCalRegion::ZSide) {
    ENUM_FOR_EACH(axis_ptr_it, FCalRegion::Axis) {
      ENUM_FOR_EACH(sign_ptr_it, FCalRegion::Sign) {
        result.insert(std::make_tuple(**side_ptr_it, **axis_ptr_it, **sign_ptr_it));
      }
    }
  }
  return result;
}

string FCalRegion::ToString(ZSide side) {
  std::map<ZSide, std::string> string_rep { {ZSide::A, "A"},
                                            {ZSide::C, "C"},
                                            {ZSide::Both, "Total"} };
  return MsgFromDict(side, string_rep);
}

string FCalRegion::ToString(Axis axis) {
  std::map<Axis, std::string> string_rep { {Axis::X, "X"},
                                           {Axis::Y, "Y"} };
  return MsgFromDict(axis, string_rep);
}

string FCalRegion::ToString(Sign sign) {
  std::map<Sign, std::string> string_rep { {Sign::Pos, "+"},
                                           {Sign::Neg, "-"} };
  return MsgFromDict(sign, string_rep);
}

/*
void FCalRegion::DumpGeoRegion(const ModuleHalf &region) {
  std::cout << ToString(std::get<0>(region))
            << ToString(std::get<1>(region))
            << ToString(std::get<2>(region)) << std::endl;
}
*/

string FCalRegion::PhiSliceFromChannel(string channel_name) {
// This is vaguely hack-y & only works with a very specific channel-naming scheme

  ZSide this_region = ZSide::A;
  if (channel_name.at(1) == '1') {
    this_region = ZSide::C;
  }

  char region_ID;
  char module_ID;
  char channel_ID;
  if (this_region == ZSide::A) {
    region_ID = 'A';
    module_ID = channel_name.at(2);
    channel_ID = channel_name.at(4);
  } else {
    region_ID = 'C';
    module_ID = channel_name.at(3);
    channel_ID = channel_name.at(5);
  }
  return PhiSliceFromID(region_ID, module_ID, channel_ID);
}

Sign FCalRegion::SignFromAxisAndPhiSlice(Axis axis, string phi_slice_name) {
  auto phi_slice_number = atoi(phi_slice_name.substr(3,2).c_str());
  if (axis == Axis::X) {
    if ((phi_slice_number < 4) || (phi_slice_number > 11)) {
      return Sign::Pos;
    } else {
      return Sign::Neg;
    }
  } else if (axis == Axis::Y) {
    if (phi_slice_number < 8) {
      return Sign::Pos;
    } else {
      return Sign::Neg;
    }
  } else {
    Error::Report("ERROR in SignFromAxisAndPhiSlice: invalid axis");
    return Sign::Pos;
  }
}
