#ifndef INCLUDE_UTIL_H_
#define INCLUDE_UTIL_H_

#include <string>

#include "error.h"
#include "void.h"

#define VPRINT(msg) if (verbose()) std::cout << msg << std::endl;

namespace Util {

Error::Expected<Void> mkdir(const std::string &path);

}

#endif
