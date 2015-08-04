#include "detector.h"

#include <string>

#include <boost/expected/expected.hpp>
#include <boost/algorithm/string.hpp>

#include "error.h"

using std::string;
using std::make_shared;

using boost::make_unexpected;
using Error::Expected;

namespace Detector {

string ToString(Name name) {
  if (name == Name::FCal) {
    return "FCal";
  }
  else {
    return "EMEC";
  }
}

// Converts a case insensitive string of "emec" or "fcal" to the appropriate
// detector enum. 
Error::Expected<Name> FromString(string&& detector_name)
{
  auto this_func_name = "DetectorFromString";
  boost::to_upper(detector_name);
  if (detector_name == "FCAL") {
    return Name::FCal;
  }
  else if (detector_name == "EMEC") {
    return Name::EMEC;
  }
  else {
    return make_unexpected(make_shared<Error::Runtime>(
        "Invalid detector name; choose \"FCal\" or \"EMEC\" (case insensitive)",
        this_func_name
    ));
  }
}

}
