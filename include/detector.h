#ifndef INCLUDE_DETECTOR_H_
#define INCLUDE_DETECTOR_H_

#include <string>

#include "error.h"

namespace Detector {

enum class Name {FCal, EMEC};

std::string ToString(Name name);
Error::Expected<Name> FromString(std::string&& detector_name);

}

#endif
