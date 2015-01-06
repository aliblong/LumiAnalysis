// With (heavy) inspiration from:
//   http://www.drdobbs.com/when-enum-just-isnt-enough-enumeration-c/184403955
#ifndef ENUM_H_
#define ENUM_H_

#include <algorithm>
#include <functional>
#include <initializer_list>
#include <iostream>
#include <set>

#define INIT_ENUM(class_name) \
  template<> \
  typename Enum<class_name>::StaticItemsSet Enum<class_name>::items_{}

#define ENUM_FOR_EACH(it_name, class_name) \
  for (auto it_name = class_name::begin(); it_name != class_name::end(); ++it_name)

template <class T>
class Enum {
 private:
  // Comparator used to determine enum item ordering
  struct ItemComparator_p {
    bool operator() (const T* e1, const T* e2) {
      return e1->value() < e2->value();
    }
  };
  struct ItemComparator_r {
    bool operator() (const T& e1, const T& e2) {
      return e1.value() < e2.value();
    }
  };

 public:
  typedef std::set<const T*, ItemComparator_p> StaticItemsSet;
  typedef typename StaticItemsSet::size_type set_size_type;
  typedef typename StaticItemsSet::const_iterator const_iterator;
  typedef typename StaticItemsSet::const_reverse_iterator const_r_iterator;
  typedef std::set<T, ItemComparator_r> Subset;

  bool operator== (const T& e2) const { return this->value() == e2.value(); }
  bool operator< (const T& e2) const { return this->value() < e2.value(); }

  int value() const { return value_; }
  static set_size_type size() { return items_.size(); }

  // Iteration
  static const_iterator begin() { return items_.cbegin(); }
  static const_iterator end() { return items_.cend(); }
  static const_r_iterator rbegin() { return items_.crbegin(); }
  static const_r_iterator rend() { return items_.crend(); }

  // Return a set which contains pointers to any number of unique enum items
  static Subset BuildSubset(std::initializer_list<T> items) {
    Subset subset;
    for (auto i: items) {
      subset.insert(i);
    }
    return subset;
  }

 protected:
  // Constructors
  explicit Enum(int value): value_(value) {
    // Warn clients if they try to add multiple items with the same value
    if (IsValidValue(value)) {
      std::cerr << "The value \'" << value << "\' has been assigned to "
                << "multiple items in this enum. All duplicates are omitted "
                << "from the set." << std::endl;
      return;
    }

    // This cast seems ugly, since the object being constructed by the
    //   inheriting class will already be of type T, but I don't see any
    //   way around this.
    items_.insert(static_cast<T*>(this));
  }
  // Compiler-generated copy and copy-assign ctors generate local copies of
  //   the static enum items

 private:
  static const T* CorrespondingEnum(int value) {
    auto EnumPredicateCorresponds = [value](const T* elem) -> bool {
      return elem->value() == value;
    };

    const_iterator it = find_if(items_.begin(),
                                items_.end(),
                                EnumPredicateCorresponds);

    return (it != items_.end()) ? *it : nullptr;
  }
  static bool IsValidValue(int value) {
    return CorrespondingEnum(value) != nullptr;
  }

  int value_;
  static StaticItemsSet items_;
};

#endif
