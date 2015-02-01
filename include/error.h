#ifndef INCLUDE_ERROR_H_
#define INCLUDE_ERROR_H_

#include <exception>
#include <iostream>
#include <memory>
#include <string>

#include <boost/expected/expected.hpp>

#include "void.h"

#define LOG_IF_ERR(result) \
  if (!result.valid()) std::cerr << result.error()->what() << std::endl;

#define RETURN_IF_ERR(result) \
  if (!result.valid()) return result.get_unexpected();

// If `result` of type Expected<T> is an error, executes `continue`.
#define CONTINUE_IF_ERR(result) \
  if (!result.valid()) { \
    LOG_IF_ERR( result ) \
    continue; \
  }

// If `result` of type Expected<T> is an error, throws a std::runtime_error.
#define THROW_IF_ERR(result)                         \
  if (!result.valid()) throw std::runtime_error(result.error()->what());

// Executes `expression` of return type Expected<T> and returns result if it
//   is an error.
#define TRY(expression)                 \
  {                                     \
    auto result_TRY = expression;           \
    RETURN_IF_ERR(result_TRY); \
  }

#define TRY_CONTINUE(expression)                 \
  {                                     \
    auto result_TRY = expression;           \
    CONTINUE_IF_ERR(result_TRY); \
  }

#define TRY_THROW(expression)                 \
  {                                     \
    auto result_TRY = expression;           \
    THROW_IF_ERR(result_TRY); \
  }

namespace Error {

class Base : public std::exception {
 public:
  ~Base() = 0;
  virtual const char* what() { return msg_.c_str(); }
 protected:
  std::string msg_;
};

template <typename T>
using Expected = boost::expected<T, std::shared_ptr<Base> >;

class System : public Base {
 public:
  System(int err_code,
         const std::string &command,
         const std::string &func_name) {
    Base::msg_ = "Error: in "+func_name+" - system command `"+command+
                   "` failed, returning error code "+std::to_string(err_code);
  }
};

class Runtime : public Base {
 public:
  Runtime(const std::string &msg,
          const std::string &func_name) {
    Base::msg_ = "Error: in "+func_name+" - "+msg;
  }
};

class File : public Base {
 public:
  File(const std::string &filepath,
       const std::string &func_name) {
    Base::msg_ = ("Error: in "+func_name+" - file `"+filepath+
                    "` was not found or could not be opened");
  }
};

class Nullptr : public Base {
 public:
  Nullptr(const std::string &func_name) {
    Base::msg_ = ("Error: "+func_name+" was passed a nullptr");
  }
};

class Logic : public Base {
 public:
  Logic(const std::string description,
        const std::string &func_name) {
    Base::msg_ = ("Error: in "+func_name+" - "+description);
  }
};

/*
// Not an extensible model since expected<T,E> doesn't play nice with pointers
enum class Type {System, File, Uninit};

class Object {
 public:
  Object () : type_(Type::Uninit), msg_("") {}
  static Object System(int err_code,
                       const std::string &command,
                       const std::string &func_name) {
    std::string msg = "Error: in "+func_name+" - system command `"+command+
                   "` failed, returning error code "+std::to_string(err_code);
    return Object(Type::System, msg);
  }
  static Object File(const std::string &filepath,
                     const std::string &func_name) {
    std::string msg = "Error: in "+func_name+" - file `"+filepath+
                         "` could not be opened";
    return Object(Type::File, msg);
  }
  const char* msg() { return msg_.c_str(); }

 private:
  Object(Type type, const std::string &msg) : type_(type), msg_(msg) {}
  Type type_;
  std::string msg_;
};
*/

}

#endif
