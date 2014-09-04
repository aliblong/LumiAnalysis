#ifndef LUMIANALYSIS_INCLUDE_JSONREADER_H_
#define LUMIANALYSIS_INCLUDE_JSONREADER_H_

#include <boost/property_tree/ptree.hpp>
#include <boost/foreach.hpp>
//#include <boost/property_tree/json_parser.hpp>
#include <string>
#include <vector>
#include <map>

class JSONReader {
  // 'using' declarations inside classes don't work with gcc apparently
  typedef typename std::string string;
  typedef typename boost::property_tree::ptree ptree;

 public:
  JSONReader(const string filename);

  template<typename T>
  T get(const string key) const {
    return pt.get<T>(key);
  }

  template<typename T>
  std::vector<T> get_vector(const string key) const {
    std::vector<T> vec;
    BOOST_FOREACH(const ptree::value_type &child, pt.get_child(key)) {
      vec.push_back(child.second.get_value<T>());
    }
    return vec;
  }

  template<typename T>
  std::map<string, T> get_map(const string key) const {
    std::map<string, T> map;
    BOOST_FOREACH(const ptree::value_type &child, pt.get_child(key)) {
      map.insert(std::make_pair(child.first, child.second.get_value<T>()));
    }
    return map;
  }

 private:
  ptree pt;
};

#endif
