#include <string>
#include <iostream>
#include <vector>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include "json_reader.h"

using std::cout;
using std::endl;
using std::string;
using std::vector;

using boost::property_tree::json_parser::read_json;

JSONReader::JSONReader(const string& filename) {
  read_json(filename, pt_);
  auto verbose = pt_.get_optional<bool>("verbose");
  if (verbose) verbose_ = *verbose;
}

void JSONReader::SetDefaults()
{
  // sick hack to read the contents of a jsonc (json with c-style comments) file into a string
  #define STRINGIFY(x) #x
  const char * jsonc_literal =
    #include "param_tree_defaults.jsonc"
    ;
  #undef STRINGIFY

  string jsonc_str(jsonc_literal);
  // erase leading and trailing parentheses (which are present to make STRINGIFY work)
  jsonc_str.erase(0, 2);
  jsonc_str.erase(jsonc_str.size()-2, 2);
  std::stringstream ss;
  ss << jsonc_str;
  read_json(ss, defaults_);
}

void JSONReader::WarnKeyNotPresent(const std::string& key) const
{
  cout << "Warning: parameter tree entry `" << key << "` missing or invalid." << endl;
}

void JSONReader::WarnDefaultNotPresent(const std::string& key) const
{
  cout << "Warning: default value for entry `" << key << "` missing or invalid." << endl;
}

/*
template<>
bool JSONReader::get<bool>(const string& key) const {
  string val = pt.get<string>(key);
  if (val == "true") {
    return true;
  }
  else if (val == "false") {
    return false;
  }
  else {
    std::cerr << "ERROR: key \'" << key << "\' corresponds to value \'"
              << val << "\', which is not a boolean. Returning false."
              << std::endl;
    return false;
  }
}
*/
