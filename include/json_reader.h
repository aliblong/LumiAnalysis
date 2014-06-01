#ifndef LUMIANALYSIS_INCLUDE_JSONREADER_H_
#define LUMIANALYSIS_INCLUDE_JSONREADER_H_

#include <boost/property_tree/ptree.hpp>
#include <boost/foreach.hpp>
//#include <boost/property_tree/json_parser.hpp>
#include <string>
#include <vector>

using std::string;
using std::vector;
using boost::property_tree::ptree;

class JSONReader {
 public:
  JSONReader(const string filename);

  template<typename T>
  T get(const string key) const {
    return pt.get<T>(key);
  }

  template<typename T>
  vector<T> get_vector(const string key) const {
    vector<T> vec;
    BOOST_FOREACH(const ptree::value_type &child, pt.get_child(key)) {
      vec.push_back(child.second.get_value<T>());
    }
    return vec;
  }

 private:
  ptree pt;
};

#endif
