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
  // Note: I couldn't figure out how to write the SFINAE part into the template parameters;
  //   turns out it ain't easy: http://rleahy.ca/blog/programming-2/variadic-sfinae/
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
      // If child_key is empty, it's an array of values (suitable for vector) rather than a subtree
      if (child_key == "") {
        auto this_val = child.second.get_value<SubTypes...>(); //SubTypes should be one type only
        result->push_back(this_val);
      }
      else {
        auto child_val = get<SubTypes...>(key+"."+child_key);
        if (!child_val) return boost::none;
        result->push_back(*child_val);
      }
    }
    return result;
  }

  /*
  template<typename FirstSubType, typename... RestOfSubTypes>
  boost::optional<typename std::enable_if<
      std::is_integral<FirstSubType>::value || std::is_floating_point<FirstSubType>::value || is_string<FirstSubType>::value,
      FirstSubType
    >::type>
    */

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

#endif
