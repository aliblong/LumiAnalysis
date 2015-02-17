#ifndef LUMIANALYSIS_INCLUDE_JSONREADER_H_
#define LUMIANALYSIS_INCLUDE_JSONREADER_H_

#include <boost/property_tree/ptree.hpp>
#include <boost/foreach.hpp>
//#include <boost/property_tree/json_parser.hpp>
#include <string>
#include <vector>
#include <map>

class JSONReader {
  typedef typename boost::property_tree::ptree ptree;

 public:
  JSONReader(const std::string& filename);

  template<typename T>
  T get(const std::string& key) const {
    return pt.get<T>(key);
  }

  template<typename T>
  std::vector<T> get_vector(const std::string& key) const {
    std::vector<T> vec;
    BOOST_FOREACH(const ptree::value_type &child, pt.get_child(key)) {
      vec.push_back(child.second.get_value<T>());
    }
    return vec;
  }

  template<typename T>
  auto get_map(const std::string& key) const {
    std::map<std::string, T> result;
    BOOST_FOREACH(const ptree::value_type &child, pt.get_child(key)) {
     result.insert(std::make_pair(child.first, child.second.get_value<T>()));
    }
    return result;
  }

 private:
  ptree pt;
};

#endif
