#ifndef LUMIANALYSIS_INCLUDE_POINT_H_
#define LUMIANALYSIS_INCLUDE_POINT_H_

#include <array>
#include <memory>
#include <vector>

#include "boost/expected/expected.hpp"

#include "error.h"

template <typename T>
using Point = std::array<T, 2>;

template <typename T>
using VectorP = std::vector<Point<T>>;

template <typename T>
Error::Expected<VectorP<T>> PointVectorFromVectors(const std::vector<T>& vec1,
                                                   const std::vector<T>& vec2) {
  // I don't want to drop in a boost header just to zip the vectors
  auto this_func_name = "PointVectorFromVectors";
  if (vec1.size() != vec2.size()) {
    auto err_msg = "vector sizes not equal";
    return boost::make_unexpected(std::make_shared<Error::Logic>(
                                    err_msg, this_func_name));
  }

  VectorP<T> result;
  auto num_points = vec1.size();
  for (size_t i = 0; i < num_points; ++i) {
    result.push_back({vec1[i], vec2[i]});
  }
  return result;
}

#endif
