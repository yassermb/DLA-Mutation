#pragma once

#include <string>
#include <vector>
#include <memory>
#include <Eigen/Sparse>

#include "cProtein.hpp"
#include "cAminoResidue.hpp"
#include "cSmoothedHistogram.hpp"


class cProteinFeaturizer {
 public:
  virtual ~cProteinFeaturizer() {}
  virtual Eigen::SparseVector<double> featurize(cProtein *protein) = 0;
};


template <class Index>
class cMultiHistogramFeaturizer : public cProteinFeaturizer {
 public:
  virtual ~cMultiHistogramFeaturizer() {}

  Eigen::SparseVector<double> featurize(cProtein *protein) final {
    accumulate_histograms(protein);
    return finalize();
  }

 protected:
  std::map<Index, std::unique_ptr<cHistogram>> features;

 private:
  virtual void accumulate_histograms(cProtein *protein) = 0;
  Eigen::SparseVector<double> finalize();
};

template <class Index>
Eigen::SparseVector<double> cMultiHistogramFeaturizer<Index>::finalize() {
  size_t dimension = 0;
  for (const auto &histogram : features) {
    dimension += histogram.second->get_vector().size();
  }

  Eigen::SparseVector<double> result(dimension);
  size_t j = 0;
  for (const auto &histogram : features) {
    for (Eigen::SparseVector<double>::InnerIterator it(histogram.second->get_vector()); it; ++it) {
      result.coeffRef(j + it.index()) = it.value();
    }
    j += histogram.second->get_vector().size();
  }
  features.clear();
  return result;
}


template <class T, class Group>
class cPairwiseHistogramFeaturizer : public cMultiHistogramFeaturizer<std::pair<Group, Group>> {
 public:
  virtual ~cPairwiseHistogramFeaturizer() {}

 protected:
  cHistogram& get_histogram(const T &first, const T &second);

 private:
  virtual void add_event(const T &first, const T &second,
                         cHistogram *histogram) = 0;
  virtual Group get_group(const T &object) const = 0;
};

template <class T, class Group>
cHistogram& cPairwiseHistogramFeaturizer<T, Group>::get_histogram(const T &first, const T &second) {
  auto it = this->features.find(std::make_pair(get_group(first), get_group(second)));
  assert(it != this->features.end());
  assert(it->second);
  return *it->second;
}


class cResiduePairwiseHist : public cPairwiseHistogramFeaturizer<cResidue,
                                        std::pair<cResidue::eType, int>> {
 public:
  // |num_DOF| - degrees of freedom
  // |num_dist_bins| - bins per degree of freedom for distance
  // |num_angle_bins| - bins per degree of freedom for angle
  // |sidechain_cutoff| - cutoff for distance between two ends of side-chains:
  //                      cutoff_CA = (R1 + R2 + sidechain_cutoff)
  // |skipped_neighbourhood| - number of neighbouring residues to skip
  // |smoothing_sigma| - sigma for truncated gaussian smoothing
  cResiduePairwiseHist(size_t num_DOF, size_t num_dist_bins, size_t num_angle_bins,
                       double sidechain_cutoff, size_t skipped_neighbourhood,
                       double smoothing_sigma);
  virtual ~cResiduePairwiseHist() {}

 private:
  double sidechain_cutoff_;
  double num_dist_bins_;
  size_t skipped_neighbourhood_;

  void add_event(const cResidue &first,
                 const cResidue &second,
                 cHistogram *histogram) final;

  typedef std::pair<cResidue::eType, int> ResidueGroup;
  ResidueGroup get_group(const cResidue &residue) const final;

  void accumulate_histograms(cProtein *protein) final;
};


class cAtomPairwiseHist : public cPairwiseHistogramFeaturizer<cAtom, std::pair<int, int>> {
 public:
  cAtomPairwiseHist(bool residue_type_dependent, bool atom_place_dependent,
                    size_t num_bins,
                    double cutoff, size_t skipped_neighbourhood,
                    double smoothing_sigma);
  virtual ~cAtomPairwiseHist() {}

 private:
  bool residue_type_dependent_;
  bool atom_place_dependent_;
  double max_interaction_dist_;
  size_t skipped_neighbourhood_;

  void add_event(const cAtom &first,
                 const cAtom &second,
                 cHistogram *histogram) final;

  void accumulate_histograms(cProtein *protein) final;
  std::pair<int, int> get_group(const cAtom &atom) const final;
  std::vector<int> get_residue_types() const;
  std::vector<int> get_atom_types() const;
};


class cHydrogenBondsHist : public cPairwiseHistogramFeaturizer<cAminoResidue, int> {
 public:
  cHydrogenBondsHist(size_t num_dist_bins, size_t num_angle_bins,
                     double max_bond_length, size_t skipped_neighbourhood,
                     double smoothing_sigma);
  virtual ~cHydrogenBondsHist() {}

 private:
  double max_interaction_dist_;
  size_t skipped_neighbourhood_;

  void add_event(const cAminoResidue &first,
                 const cAminoResidue &second,
                 cHistogram *histogram) final;

  int get_group(const cAminoResidue &atom) const final { return 0; }
  void accumulate_histograms(cProtein *protein) final;
};


class cSolvationShellHist : public cPairwiseHistogramFeaturizer<cResidue, int> {
 public:
  cSolvationShellHist(double sa_radius,
                      double sa_smoothed_margin,
                      size_t num_dist_bins,
                      size_t num_angle_bins,
                      double cutoff,
                      double smoothing_sigma);
  virtual ~cSolvationShellHist() {}

 private:
  double sa_radius_;
  double sa_smoothed_margin_;
  double max_interaction_dist_;

  void add_event(const cResidue &residue,
                 const cResidue &water,
                 cHistogram *histogram) final;

  int get_group(const cResidue &residue) const final;
  void accumulate_histograms(cProtein *protein) final;
};
