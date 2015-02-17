#ifndef LUMIANALYSIS_INCLUDE_ATLAS_LABEL_OPTIONS_H_
#define LUMIANALYSIS_INCLUDE_ATLAS_LABEL_OPTIONS_H_

class ATLASLabelOptions {
 public:
  ATLASLabelOptions(std::string subheading) : subheading_(std::move(subheading)) {}
  // TODO: write this implementation
  //ATLASLabelOptions(std::string params_filepath);
  ~ATLASLabelOptions(){};

  auto& text_size(Float_t text_size) { text_size_ = text_size; return *this; }

  auto& x(Float_t x) { x_ = x; return *this; }
  auto& y(Float_t y) { y_ = y; return *this; }
  auto& subheading_offset(Float_t subheading_offset) {
    subheading_offset_ = subheading_offset; return *this;
  }

  const auto& subheading() const { return subheading_; }

  auto text_size() const { return text_size_; }

  auto use_NDC() const { return use_NDC_; }
  auto x() const { return x_; }
  auto y() const { return y_; }
  auto subheading_offset() const { return subheading_offset_; }

 private:
  std::string subheading_;

  Float_t text_size_ = 0.04;

  bool use_NDC_ = true;
  Float_t x_ = 0.65;
  Float_t y_ = 0.84;
  Float_t subheading_offset_ = 0.10;
};

#endif
