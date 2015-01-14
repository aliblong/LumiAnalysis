#ifndef INCLUDE_UTIL_H_
#define INCLUDE_UTIL_H_

#include <string>

#include "error.h"
#include "void.h"

namespace Util {

Error::Expected<Void> mkdir(const std::string &path);

}

#endif
