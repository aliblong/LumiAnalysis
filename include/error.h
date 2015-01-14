#ifndef INCLUDE_ERROR_H_
#define INCLUDE_ERROR_H_

#include <exception>
#include <memory>
#include <string>

#include <boost/expected/expected.hpp>

#define TRY(expression) \
  if (!expression.valid()) return expression

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

class FileNotOpen : public Base {
 public:
  FileNotOpen(const std::string &filepath,
              const std::string &func_name) {
    Base::msg_ = ("Error: in "+func_name+" - file `"+filepath+
                    "` could not be opened");
  }
};

class Nullptr : public Base {
 public:
  Nullptr(const std::string &func_name) {
    Base::msg_ = ("Error: "+func_name+" was passed a nullptr");
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

void Report(const std::string &err_msg);

void Report(const std::string &err_msg, unsigned indent_level);

}

#endif
