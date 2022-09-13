#include "cSmoothedHistogram.hpp"

#include <iostream>
#include <functional>
#include <numeric>


const double kWindowRatio = 1.0 / 2;
const double kEps = 1e-9;


cSmoothedHistogram::cSmoothedHistogram(const std::vector<Dimension> &dimensions)
    : dimensions_(dimensions) {
  assert(dimensions_.size());
  size_t total_size = std::accumulate(
    dimensions_.begin(), dimensions_.end(), 1,
    [](size_t size, const Dimension &dim) { return size * dim.size(); }
  );
  vector_ = Eigen::SparseVector<double>(total_size);
}


double CenteredGaussianCDF(double x, double sigma) {
  return (1 + std::erf(x / (sigma * std::sqrt(2.0))));
}

double TruncatedGaussianCDF(double x,
                            double left_border, double right_border,
                            double sigma = 1.0) {
  if (x <= left_border)
    return 0.0;
  if (x >= right_border)
    return 1.0;

  double a = CenteredGaussianCDF(left_border, sigma);
  double b = CenteredGaussianCDF(right_border, sigma);
  return (CenteredGaussianCDF(x, sigma) - a) / (b - a);
}

double IntegrateTruncatedCenteredGaussian(double from, double to,
                                          double sigma, double support) {
  // To avoid floating point pitfalls and losing of sparsity
  if (sigma < kEps) {
    return from <= 0 && to > 0;
  }

  return TruncatedGaussianCDF(to, -support / 2, support / 2, sigma)
          - TruncatedGaussianCDF(from, -support / 2, support / 2, sigma);
}


cSmoothedHistogram::Dimension::Dimension(double left_border, eBorderType left_border_type,
                                         double right_border, eBorderType right_border_type,
                                         size_t size, double smoothing_sigma)
      : left_border_(left_border), left_border_type_(left_border_type),
        right_border_(right_border), right_border_type_(right_border_type),
        smoothing_sigma_(smoothing_sigma) {
  assert(left_border < right_border);
  if (left_border_type_ == CYCLIC || right_border_type_ == CYCLIC) {
    assert(left_border_type_ == right_border_type_);
  }

  bins_.reserve(size);
  for (size_t i = 0; i < size; ++i) {
    bins_.push_back({ static_cast<double>(i) / size,
                      static_cast<double>(i + 1) / size });
  }
}

std::vector<double> cSmoothedHistogram::Dimension::Integrate(double x) const {
  x = (x - left_border_) / (right_border_ - left_border_);

  std::vector<double> integrals(bins_.size(), 0);

  double support = (bins_[0].end - bins_[0].begin) * kWindowRatio;

  // (..., [0, 1, 2, 3, 4, 5], ...)
  for (size_t i = 0; i < bins_.size(); ++i) {
    double integral = IntegrateTruncatedCenteredGaussian(
      bins_[i].begin - x,
      bins_[i].end - x,
      smoothing_sigma_ / (right_border_ - left_border_),
      support
    );
    integrals[i] += integral;
  }

  double hist_length = bins_.back().end - bins_.front().begin;

  // (0, 1, 2, 3, 4, 5, [0, 1, 2, 3, 4, 5], 0, 1, 2, 3, 4, 5)
  if (left_border_type_ == CYCLIC) {
    for (size_t i = 0; i < bins_.size(); ++i) {
      double integral = IntegrateTruncatedCenteredGaussian(
        bins_[i].begin + hist_length - x,
        bins_[i].end + hist_length - x,
        smoothing_sigma_,
        support
      );
      integrals[i] += integral;
    }
    for (size_t i = 0; i < bins_.size(); ++i) {
      double integral = IntegrateTruncatedCenteredGaussian(
        bins_[i].begin - hist_length - x,
        bins_[i].end - hist_length - x,
        smoothing_sigma_,
        support
      );
      integrals[i] += integral;
    }
  }

  // (..., 5, 4, 3, 2, 1, 0, [0, 1, 2, 3, 4, 5, ...)
  if (left_border_type_ == CLOSED) {
    for (size_t i = 0; i < bins_.size(); ++i) {
      double integral = IntegrateTruncatedCenteredGaussian(
        bins_.front().begin - (bins_[i].end - bins_.front().begin) - x,
        bins_.front().begin - (bins_[i].begin - bins_.front().begin) - x,
        smoothing_sigma_,
        support
      );
      integrals[i] += integral;
    }
  }
  // (..., 0, 1, 2, 3, 4, 5], 5, 4, 3, 2, 1, 0, ...)
  if (right_border_type_ == CLOSED) {
    for (size_t i = 0; i < bins_.size(); ++i) {
      double integral = IntegrateTruncatedCenteredGaussian(
        bins_.back().end + hist_length - (bins_[i].end - bins_.front().begin) - x,
        bins_.back().end + hist_length - (bins_[i].begin - bins_.front().begin) - x,
        smoothing_sigma_,
        support
      );
      integrals[i] += integral;
    }
  }
  return integrals;
}

size_t cSmoothedHistogram::convert_tensor_idx(const std::vector<size_t> &indices) const {
  assert(indices.size() == dimensions_.size());

  size_t idx = 0;
  for (size_t i = 0; i < dimensions_.size(); ++i) {
    assert(indices[i] < dimensions_[i].size());
    idx = idx * dimensions_[i].size() + indices[i];
  }
  return idx;
}


template <typename T>
std::vector<std::vector<T>> CartesianProduct(const std::vector<std::vector<T>> &mesh) {
  if (std::any_of(mesh.begin(), mesh.end(),
                  [](const std::vector<T> &vector) { return vector.size() == 0; })) {
    return std::vector<std::vector<T>>();
  }

  size_t dimension = mesh.size();
  size_t num_points = 1;
  for (const auto &vector : mesh) {
    num_points *= vector.size();
  }
  std::vector<std::vector<T>> result(num_points, std::vector<T>(dimension));

  for (size_t i = 0, mod = 1; i < dimension; ++i) {
    for (size_t j = 0; j < result.size(); ++j) {
      result[j][i] = mesh[i][(j / mod) % mesh[i].size()];
    }
    mod *= mesh[i].size();
  }
  return result;
};

void cSmoothedHistogram::add_event(const std::vector<double> &point, double weight) {
  assert(point.size());
  assert(point.size() == dimensions_.size());

  std::vector<std::vector<size_t>> indices_mesh(dimensions_.size());
  std::vector<std::vector<double>> integrals_mesh(dimensions_.size());

  for (size_t i = 0; i < dimensions_.size(); ++i) {
    auto integrals = dimensions_[i].Integrate(point[i]);
    assert(integrals.size() == dimensions_[i].size());

    for (size_t bin = 0; bin < integrals.size(); ++bin) {
      if (integrals[bin] > 0.0) {
        indices_mesh[i].push_back(bin);
        integrals_mesh[i].push_back(integrals[bin]);
      }
    }
  }
  auto indices = CartesianProduct(indices_mesh);
  auto integrals = CartesianProduct(integrals_mesh);

  assert(indices.size() == integrals.size());

  for (size_t i = 0; i < indices.size(); ++i) {
    double integral = std::accumulate(integrals[i].begin(), integrals[i].end(),
                                       1.0, std::multiplies<double>());
    if (integral > 0.0) {
      vector_.coeffRef(convert_tensor_idx(indices[i])) += integral * weight;
    }
  }
}
