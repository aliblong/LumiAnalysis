#ifndef LUMIANALYSIS_INCLUDE_JSONREADER_H_
#define LUMIANALYSIS_INCLUDE_JSONREADER_H_

#include <boost/property_tree/ptree.hpp>
#include <boost/foreach.hpp>
//#include <boost/property_tree/json_parser.hpp>
#include <string>
#include <vector>
#include <iostream>

#include "boost/container/flat_map.hpp"
#include "boost/optional.hpp"

class JSONReader {
  typedef typename boost::property_tree::ptree ptree;

 public:
  JSONReader(const std::string& filename);

  template<typename T>
  auto get(const std::string& key) const {
    return pt.get_optional<T>(key);
  }

  template<typename T>
  T get(const std::string& key, T default_val) const {
    return pt.get<T>(key, default_val);
  }

  template<typename T>
  T get(const std::string& key, T default_val, bool verbose) const {
    auto result = pt.get_optional<T>(key);
    if (!result) {
      if (verbose) {
        std::cout << "Warning: parameter tree entry `" << key << 
          "` missing or invalid. Defaulting to value `" << default_val << "`." <<std::endl;
      }
      return default_val;
    }
    return *result;
  }

  template<typename T>
  auto get_vector(const std::string& key) const {
    boost::optional<std::vector<T>> result;
    auto node = pt.get_child_optional(key);
    if (!node) return result;
    result = std::vector<T>();
    BOOST_FOREACH(const ptree::value_type &child, *node) {
      result->push_back(child.second.get_value<T>());
    }
    return result;
  }

  template<typename T>
  auto get_vector(const std::string& key, std::vector<T>&& default_val) const {
    auto node = pt.get_child_optional(key);
    if (!node) return default_val;
    std::vector<T> result;
    BOOST_FOREACH(const ptree::value_type &child, *node) {
      result.push_back(child.second.get_value<T>());
    }
    return result;
  }

  template<typename T>
  auto get_vector(const std::string& key, std::vector<T>&& default_val, bool verbose) const {
    auto node = pt.get_child_optional(key);
    if (!node) {
      if (verbose) {
        std::cout << "Warning: parameter tree entry `" << key << 
          "` missing or invalid." <<std::endl;
      }
      return default_val;
    }
    std::vector<T> result;
    BOOST_FOREACH(const ptree::value_type &child, *node) {
      result.push_back(child.second.get_value<T>());
    }
    return result;
  }

  template<typename T>
  auto get_map(const std::string& key) const {
    boost::optional<boost::container::flat_map<std::string, T>> result;
    auto node = pt.get_child_optional(key);
    if (!node) return result;
    result = boost::container::flat_map<std::string, T>();
    BOOST_FOREACH(const ptree::value_type &child, *node) {
     result->insert(std::make_pair(child.first, child.second.get_value<T>()));
    }
    return result;
  }

  template<typename T>
  auto get_map(
      const std::string& key,
      boost::container::flat_map<std::string,T>&& default_val
  ) const {
    auto node = pt.get_child_optional(key);
    if (!node) return default_val;
    boost::container::flat_map<std::string, T> result;
    BOOST_FOREACH(const ptree::value_type &child, *node) {
     result.insert(std::make_pair(child.first, child.second.get_value<T>()));
    }
    return result;
  }

  template<typename T>
  auto get_map(
      const std::string& key,
      boost::container::flat_map<std::string,T>&& default_val,
      bool verbose
  ) const {
    auto node = pt.get_child_optional(key);
    if (!node) {
      if (verbose) {
        std::cout << "Warning: parameter tree entry `" << key << 
          "` missing or invalid." <<std::endl;
      }
      return default_val;
    }
    boost::container::flat_map<std::string, T> result;
    BOOST_FOREACH(const ptree::value_type &child, *node) {
     result.insert(std::make_pair(child.first, child.second.get_value<T>()));
    }
    return result;
  }

  template<typename T>
  auto get_map_of_vectors(const std::string& key) const {
    boost::optional<boost::container::flat_map<std::string, std::vector<T>>> result;
    auto node = pt.get_child_optional(key);
    if (!node) return result;
    result = boost::container::flat_map<std::string, T>();
    for (const auto& child: pt.get_child(key)) {
      auto key_to_insert = child.first;
      auto val_to_insert = get_vector<T>(key+key_to_insert);
      result->insert(std::make_pair(key_to_insert, val_to_insert));
    }
    return result;
  }

  template<typename T>
  auto get_map_of_vectors(
      const std::string& key,
      boost::container::flat_map<std::string, std::vector<T>>&& default_val
  ) const {
    boost::container::flat_map<std::string, std::vector<T>> result;
    auto node = pt.get_child_optional(key);
    if (!node) return default_val;
    for (const auto& child: *node) {
      auto key_to_insert = child.first;
      auto val_to_insert = get_vector<T>(key+key_to_insert);
      result.insert(std::make_pair(key_to_insert, val_to_insert));
    }
    return result;
  }

  template<typename T>
  auto get_map_of_vectors(
      const std::string& key,
      boost::container::flat_map<std::string, std::vector<T>>&& default_val,
      bool verbose
  ) const {
    boost::container::flat_map<std::string, std::vector<T>> result;
    auto node = pt.get_child_optional(key);
    if (!node) {
      if (verbose) {
        std::cout << "Warning: parameter tree entry `" << key << 
          "` missing or invalid." <<std::endl;
      }
      return default_val;
    }
    for (const auto& child: *node) {
      auto key_to_insert = child.first;
      auto val_to_insert = get_vector<T>(key+key_to_insert, {}, true);
      result.insert(std::make_pair(key_to_insert, val_to_insert));
    }
    return result;
  }

 private:
  ptree pt;
};

#endif
