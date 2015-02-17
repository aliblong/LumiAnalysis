#include "fcal_region.h"

#include <iostream>
#include <map>
#include <sstream>
#include <string>

#include "boost/expected/expected.hpp"

#include "error.h"

using FCalRegion::ZSide;
using FCalRegion::Axis;
using FCalRegion::Sign;
using FCalRegion::ModuleHalfSet;

using std::string;
using std::make_shared;

using boost::make_unexpected;
using Error::Expected;

namespace {

struct LineID {
  LineID(char r, char m, char c)
    : region_ID(r),
      module_ID(m),
      channel_ID(c) {}
  char region_ID;
  char module_ID;
  char channel_ID;
};

const std::map<ZSide, std::string> Z_SIDE_STR_REP { {ZSide::A, "A"},
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
  }
  else {
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

Expected<string> PhiSliceFromID(const LineID line_ID) {
  auto this_func_name = "PhiSliceFromID";

  auto module_ID_int = IntFromChar(line_ID.module_ID);
  if ( (module_ID_int < 0) || (module_ID_int > 9) ) {
    auto err_msg = "module_ID input falls outside [0,9]";
    return make_unexpected(make_shared<Error::Logic>(err_msg, this_func_name));
  }
  auto channel_ID_int = IntFromChar(line_ID.channel_ID);
  if ( (channel_ID_int < 0) || (channel_ID_int > 9) ) {
    auto err_msg = "channel_ID input falls outside [0,9]";
    return make_unexpected(make_shared<Error::Logic>(err_msg, this_func_name));
  }

  unsigned phi_slice_ID_int = (module_ID_int * 2) + (channel_ID_int / 4);
  string region_ID_str(1, line_ID.region_ID);
  string phi_slice = region_ID_str + "1." + StringFromInt(phi_slice_ID_int, 2);
  return phi_slice;
}

}

ModuleHalfSet FCalRegion::CreateModuleHalfSet() {
  ModuleHalfSet result;
  for (auto side: FCalRegion::Z_SIDES) {
    for (auto axis: FCalRegion::AXES) {
      for (auto sign: FCalRegion::SIGNS) {
        result.insert(std::make_tuple(side, axis, sign));
      }
    }
  }

  return result;
}

string FCalRegion::ToString(ZSide side) {
  return MsgFromDict(side, Z_SIDE_STR_REP);
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

Expected<string> FCalRegion::PhiSliceFromChannel(string channel_name) {
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
  }
  else {
    region_ID = 'C';
    module_ID = channel_name.at(3);
    channel_ID = channel_name.at(5);
  }
  return PhiSliceFromID(LineID(region_ID, module_ID, channel_ID));
}

Sign FCalRegion::SignFromAxisAndPhiSlice(Axis axis, string phi_slice_name) {
  auto phi_slice_number = atoi(phi_slice_name.substr(3,2).c_str());
  if (axis == Axis::X) {
    if ((phi_slice_number < 4) || (phi_slice_number > 11)) {
      return Sign::Pos;
    }
    else {
      return Sign::Neg;
    }
  }
  else { //if (axis == Axis::Y) {
    if (phi_slice_number < 8) {
      return Sign::Pos;
    }
    else {
      return Sign::Neg;
    }
  }
}

// A half-assed test
bool FCalRegion::IsValidChannel(const string& channel_name) {
  auto size = channel_name.size();
  if (size < 5 || size > 6) return false;
  if (channel_name[0] != 'M') return false;
  if (channel_name[1] != '1' && channel_name[1] != '8') return false;
  return true;
}
