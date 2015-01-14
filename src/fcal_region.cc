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
using FCalRegion::ModuleHalfSet;

namespace {

const std::map<ZSide, std::string> z_side_str_rep { {ZSide::A, "A"},
                                                    {ZSide::C, "C"},
                                                    {ZSide::Both, "Total"} };
const std::map<Axis, std::string> axis_str_rep { {Axis::X, "X"},
                                                 {Axis::Y, "Y"} };
const std::map<Sign, std::string> sign_str_rep { {Sign::Pos, "+"},
                                                 {Sign::Neg, "-"} };

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

ModuleHalfSet FCalRegion::CreateModuleHalfSet() {
  ModuleHalfSet result;
  for (const auto &side: FCalRegion::Z_SIDES) {
    for (const auto &axis: FCalRegion::AXES) {
      for (const auto &sign: FCalRegion::SIGNS) {
        result.insert(std::make_tuple(side, axis, sign));
      }
    }
  }

  return result;
}

string FCalRegion::ToString(ZSide side) {
  return MsgFromDict(side, z_side_str_rep);
}

string FCalRegion::ToString(Axis axis) {
  return MsgFromDict(axis, axis_str_rep);
}

string FCalRegion::ToString(Sign sign) {
  return MsgFromDict(sign, sign_str_rep);
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
