#include <string>

namespace Error {

void Report(const std::string &err_msg);

void Report(const std::string &err_msg, unsigned indent_level);

}
