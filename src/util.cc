#include <memory>
#include <string>

#include "boost/expected/expected.hpp"

#include "error.h"
#include "util.h"
#include "void.h"

using std::string;

Error::Expected<Void> Util::mkdir(const string &path) {
  const char *this_func_name = "Util::mkdir";
  const string mkdir_command = "mkdir -p " + path;

  int err = system( mkdir_command.c_str() );
  if (err) {
    return boost::make_unexpected(std::make_shared<Error::System>(err, mkdir_command, this_func_name));
  }

  return Void();
}
