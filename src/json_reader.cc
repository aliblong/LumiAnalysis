#include <string>
#include <iostream>
#include <vector>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include "json_reader.h"

using std::string;
using std::vector;

using boost::property_tree::ptree;

JSONReader::JSONReader(const string filename) {
  boost::property_tree::json_parser::read_json(filename, pt);
}

template<>
bool JSONReader::get<bool>(const string key) const {
  string val = pt.get<string>(key);
  if (val == "true") {
    return true;
  } else if (val == "false") {
    return false;
  } else {
    std::cerr << "ERROR: key \'" << key << "\' corresponds to value \'"
              << val << "\', which is not a boolean. Returning false."
              << std::endl;
    return false;
  }
}
