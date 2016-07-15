#ifndef LUMIANALYSIS_INCLUDE_JSONREADER_H_
#define LUMIANALYSIS_INCLUDE_JSONREADER_H_

#include <boost/property_tree/ptree.hpp>
#include <boost/foreach.hpp>
//#include <boost/property_tree/json_parser.hpp>
#include <string>
#include <vector>

#include "boost/container/flat_map.hpp"

class JSONReader {
  typedef typename boost::property_tree::ptree ptree;

 public:
  JSONReader(const std::string& filename);

  template<typename T>
  T get(const std::string& key) const {
    return pt.get<T>(key);
  }

  template<typename T>
  T get(const std::string& key, T default_val) const {
    return pt.get<T>(key, default_val);
  }

  template<typename T>
  auto get_vector(const std::string& key) const {
    std::vector<T> vec;
    BOOST_FOREACH(const ptree::value_type &child, pt.get_child(key)) {
      vec.push_back(child.second.get_value<T>());
    }
    return vec;
  }

  template<typename T>
  auto get_vector(const std::string& key, std::vector<T>&& default_val) const {
    std::vector<T> vec;
    auto node = pt.get_child_optional(key);
    if (!node) {
      return default_val;
    }
    else {
      BOOST_FOREACH(const ptree::value_type &child, *node) {
        vec.push_back(child.second.get_value<T>());
      }
    return vec;
    }
  }

  template<typename T>
  auto get_map(const std::string& key) const {
    boost::container::flat_map<std::string, T> result;
    BOOST_FOREACH(const ptree::value_type &child, pt.get_child(key)) {
     result.insert(std::make_pair(child.first, child.second.get_value<T>()));
    }
    return result;
  }

  template<typename T>
  auto get_map(const std::string& key, boost::container::flat_map<std::string,T>&& default_val) const {
    boost::container::flat_map<std::string, T> result;
    auto node = pt.get_child_optional(key);
    if (!node) {
      return default_val;
    }
    else {
      BOOST_FOREACH(const ptree::value_type &child, *node) {
       result.insert(std::make_pair(child.first, child.second.get_value<T>()));
      }
      return result;
    }
  }

  template<typename T>
  auto get_map_of_vectors(const std::string& key) const {
    boost::container::flat_map<std::string, std::vector<T>> result;
    for (const auto& child: pt.get_child(key)) {
      auto key_to_insert = child.first;
      auto val_to_insert = get_vector<T>(key+key_to_insert);
      result.insert(std::make_pair(key_to_insert, val_to_insert));
    }
    return result;
  }

  template<typename T>
  auto get_map_of_vectors(const std::string& key, boost::container::flat_map<std::string, std::vector<T>>&& default_val) const {
    boost::container::flat_map<std::string, std::vector<T>> result;
    auto node = pt.get_child_optional(key);
    if (!node) {
      return default_val;
    }
    else {
      for (const auto& child: *node) {
        auto key_to_insert = child.first;
        auto val_to_insert = get_vector<T>(key+key_to_insert);
        result.insert(std::make_pair(key_to_insert, val_to_insert));
      }
    }
    return result;
  }

 private:
  ptree pt;
};

#endif
