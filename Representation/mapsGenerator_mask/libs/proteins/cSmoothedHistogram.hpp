#pragma once

#include <vector>
#include <Eigen/Sparse>


class cHistogram {
 public:
  virtual ~cHistogram() {}

  virtual void add_event(const std::vector<double> &point, double weight = 1) = 0;
  virtual size_t num_dims() const = 0;

  virtual Eigen::SparseVector<double>& get_vector() = 0;
  virtual const Eigen::SparseVector<double>& get_vector() const = 0;
};


class cDummyHistogram : public cHistogram {
 public:
  virtual ~cDummyHistogram() {}

  virtual void add_event(const std::vector<double> &point, double weight = 1) {}
  virtual size_t num_dims() const { return 0; };

  virtual Eigen::SparseVector<double>& get_vector() { return vector_; }
  virtual const Eigen::SparseVector<double>& get_vector() const { return vector_; }

 protected:
  Eigen::SparseVector<double> vector_;
};


// Histogram, where events are smoothed by separable smoothing kernel
// So, filter is applied to each dimension separately
class cSmoothedHistogram : public cDummyHistogram {
 public:
  class Dimension;

  explicit cSmoothedHistogram(const std::vector<Dimension> &dimensions);
  virtual ~cSmoothedHistogram() {}

  virtual void add_event(const std::vector<double> &point, double weight = 1) final;
  virtual size_t num_dims() const { return dimensions_.size(); };

 private:
  std::vector<Dimension> dimensions_;

  // (i1, i2, i3) -> i3 + d3*(i2 + d2*i1)
  size_t convert_tensor_idx(const std::vector<size_t> &point) const;
};


class cSmoothedHistogram::Dimension {
 public:
  // The border type of a dimension assigns smoothing method
  // and it is determined by physical interpretation of event
  enum eBorderType {
    OPEN,  // the event can happen beyond the border
    CLOSED,  // the event can not happen beyond the border
    CYCLIC  // Both dimension's borders coincide
  };

  Dimension(double left_border, eBorderType left_border_type,
            double right_border, eBorderType right_border_type,
            size_t size, double smoothing_sigma);

  std::vector<double> Integrate(double x) const;

  size_t size() const { return bins_.size(); }

 private:
  double left_border_;
  eBorderType left_border_type_;
  double right_border_;
  eBorderType right_border_type_;

  double smoothing_sigma_;

  struct Bin {
    double begin;
    double end;
  };
  std::vector<Bin> bins_;
};
