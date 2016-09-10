#ifndef LUMIANALYSIS_INCLUDE_EXTENDED_TYPE_TRAITS_H_
#define LUMIANALYSIS_INCLUDE_EXTENDED_TYPE_TRAITS_H_

#include <string>
#include <vector>

#include "boost/container/flat_map.hpp"

template<class T>
struct is_string {
    static bool const value = false;
};

template<>
struct is_string<std::string> {
    static bool const value = true;
};

template<>
struct is_string<const char*> {
    static bool const value = true;
};

template<class T>
struct is_vector {
    static bool const value = false;
};

template<class T>
struct is_vector<std::vector<T>> {
    static bool const value = true;
};

template<class T>
struct is_map {
    static bool const value = false;
};

template<class K, class V>
struct is_map<boost::container::flat_map<K, V>> {
    static bool const value = true;
};

#endif
