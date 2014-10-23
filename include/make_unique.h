#ifndef LUMIANALYSIS_INCLUDE_MAKE_UNIQUE_H_
#define LUMIANALYSIS_INCLUDE_MAKE_UNIQUE_H_

#include <memory>
#include <utility>

// Implementation of std::make_unique since we don't have access to C++14
namespace std {

  template<typename T, typename... Args>
  std::unique_ptr<T> make_unique(Args&&... args) {
      return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
  }

}

#endif
