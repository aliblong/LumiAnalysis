#include "error.h"

#include <iostream>
#include <string>

using std::cerr;
using std::endl;

using std::string;

void Error::Report(const string &err_msg) {
  Error::Report(err_msg, 0);
}

void Error::Report(const string &err_msg, unsigned indent_level) {
  const unsigned max_indent = 9;
  if (indent_level > max_indent) {
    cerr << "WARNING: too many (>" << max_indent
         << ") indent levels in your error message." << endl;
    indent_level = max_indent;
  }
  string out_msg;
  for (int i = 0; i < indent_level; ++i) out_msg += '\t';
  cerr << err_msg << endl;
}
