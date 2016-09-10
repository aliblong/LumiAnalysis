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

#include "extended_type_traits.h"

class JSONReader {
  typedef typename boost::property_tree::ptree ptree;

 public:
  JSONReader(const std::string& filename);
  void SetDefaults();
  bool verbose() const { return verbose_; }

  template<typename T>
  boost::optional<typename std::enable_if<
      std::is_integral<T>::value || std::is_floating_point<T>::value || is_string<T>::value,
      T
    >::type>
  get(const std::string& key) const {
    auto result = pt_.get_optional<T>(key);
    if (!result && !defaults_.empty()) {
      result = defaults_.get_optional<T>(key);
      if (verbose_) {
        WarnKeyNotPresent(key);
        if (result) {
          WarnDefaultUsed(*result);
        }
        else {
          WarnDefaultNotPresent(key);
        }
      }
    }
    return result;
  }

  // T is the entire vector type; Subtypes are the subsequent nested types
  // e.g. for T = vector<vector<int>>, Subtypes are vector<int> and int (in that order!)
  template<typename T, typename... SubTypes>
  boost::optional<typename std::enable_if<is_vector<T>::value, T>::type>
  get(const std::string& key) const {
    boost::optional<T> result;
    auto node = pt_.get_child_optional(key);
    if (!node && !defaults_.empty()) node = defaults_.get_child_optional(key);
    if (!node) return result;
    result = T();
    BOOST_FOREACH(const ptree::value_type &child, *node) {
      auto child_key = child.first;
      auto child_val = get<SubTypes...>(key+"."+child_key);
      if (!child_val) return boost::none;
      result->push_back(*child_val);
    }
    return result;
  }

  // T is the entire vector type; Subtypes are the subsequent nested types, noting that key type
  //   is always string
  // e.g. for T = map<string, vector<float>>, Subtypes are vector<float> and float (in that order!)
  template<typename T, typename... SubTypes>
  boost::optional<typename std::enable_if<is_map<T>::value, T>::type>
  get(const std::string& key) const {
    boost::optional<T> result;
    auto node = pt_.get_child_optional(key);
    if (!node) return result;
    result = T();
    BOOST_FOREACH(const ptree::value_type &child, *node) {
      auto child_key = child.first;
      auto child_val = get<SubTypes...>(key+"."+child_key);
      if (!child_val) return boost::none;
     result->insert(std::make_pair(child_key, *child_val));
    }
    return result;
  }

 private:
  void WarnKeyNotPresent(const std::string& key) const;
  template<typename T>
  void WarnDefaultUsed(T default_val) const {
    std::cout << "Warning: defaulting to value `" << default_val << "`." << std::endl;
  }
  void WarnDefaultNotPresent(const std::string& key) const;

  ptree pt_;
  ptree defaults_;
  bool verbose_ = true;
};

/*
  template<typename T>
  auto get(const std::string& key) const {
    auto result = pt_.get_optional<T>(key);
    if (!result) result = defaults_.get_optional<T>(key);
    return result;
  }

  template<typename T>
  T get(const std::string& key, T default_val) const {
    return pt_.get<T>(key, default_val);
  }

  template<typename T>
  T get(const std::string& key, T default_val, bool verbose) const {
    auto result = pt_.get_optional<T>(key);
    if (!result) {
      if (verbose) {
      }
      return default_val;
    }
    return *result;
  }
  template<typename T>
  auto get_vector(const std::string& key, std::vector<T>&& default_val) const {
    auto node = pt_.get_child_optional(key);
    if (!node) return default_val;
    std::vector<T> result;
    BOOST_FOREACH(const ptree::value_type &child, *node) {
      result.push_back(child.second.get_value<T>());
    }
    return result;
  }

  template<typename T>
  auto get_vector(const std::string& key, std::vector<T>&& default_val, bool verbose) const {
    auto node = pt_.get_child_optional(key);
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
  auto get_map_of_vectors(
      const std::string& key,
      boost::container::flat_map<std::string, std::vector<T>>&& default_val
  ) const {
    boost::container::flat_map<std::string, std::vector<T>> result;
    auto node = pt_.get_child_optional(key);
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
    auto node = pt_.get_child_optional(key);
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

  template<typename T>
  auto get_map(
      const std::string& key,
      boost::container::flat_map<std::string,T>&& default_val
  ) const {
    auto node = pt_.get_child_optional(key);
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
    auto node = pt_.get_child_optional(key);
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
    auto node = pt_.get_child_optional(key);
    if (!node) return result;
    result = boost::container::flat_map<std::string, T>();
    for (const auto& child: pt_.get_child(key)) {
      auto key_to_insert = child.first;
      auto val_to_insert = get_vector<T>(key+key_to_insert);
      result->insert(std::make_pair(key_to_insert, val_to_insert));
    }
    return result;
  }


*/

#endif
