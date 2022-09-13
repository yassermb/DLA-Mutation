#include "cProteinFeaturizer.hpp"

#include <assert.h>
#include <map>

#include "cAminoResidue.hpp"
#include "mathFunctions.hpp"
#include "cGrid.hpp"
#include "cSmoothedHistogram.hpp"
#include "NeighbourSearch.hpp"
#include "SolventGrid.hpp"



cResiduePairwiseHist::cResiduePairwiseHist(size_t num_DOF,
                                           size_t num_dist_bins, size_t num_angle_bins,
                                           double sidechain_cutoff, size_t skipped_neighbourhood,
                                           double smoothing_sigma)
    : sidechain_cutoff_(sidechain_cutoff),
      num_dist_bins_(num_dist_bins),
      skipped_neighbourhood_(skipped_neighbourhood) {
  std::vector<ResidueGroup> groups = {
    ResidueGroup(cResidue::NUCLEIC, 0),
    ResidueGroup(cResidue::WATER, 0),
    ResidueGroup(cResidue::ION, 0),
  };

  for (int amino_type = cAminoResidue::ALA;
            amino_type < cAminoResidue::AMINO_TYPE_END; ++amino_type) {
    groups.push_back(ResidueGroup(cResidue::AMINO, amino_type));
  }

  for (const auto &first_group : groups) {
    for (const auto &second_group : groups) {
      if (first_group > second_group)
        continue;

      if (first_group.first != cResidue::AMINO
          || second_group.first != cResidue::AMINO) {
        features[std::make_pair(first_group, second_group)] = std::unique_ptr<cHistogram>(new cDummyHistogram());
        continue;
      }

      double residue_interaction_cutoff = std::max(1e-20,
        contact_dist(static_cast<cAminoResidue::eAminoType>(first_group.second),
                     static_cast<cAminoResidue::eAminoType>(second_group.second))
        + sidechain_cutoff_);

      double displacement_smoothing_sigma = smoothing_sigma * std::sqrt(2);
      double dihedral_smoothing_sigma = smoothing_sigma * std::sqrt(2);
      double theta_smoothing_sigma = smoothing_sigma;
      // TODO: dependency on theta
      double phi_smoothing_sigma = smoothing_sigma;

      typedef cSmoothedHistogram::Dimension Dim;
      switch (num_DOF) {
        case 1:
          features[std::make_pair(first_group, second_group)]
                = std::unique_ptr<cHistogram>(new cSmoothedHistogram({
            { 0, Dim::CLOSED, residue_interaction_cutoff, Dim::OPEN, num_dist_bins, displacement_smoothing_sigma }
          }));
          break;
        case 2:
          features[std::make_pair(first_group, second_group)]
                = std::unique_ptr<cHistogram>(new cSmoothedHistogram({
            { 0, Dim::CLOSED, residue_interaction_cutoff, Dim::OPEN, num_dist_bins, displacement_smoothing_sigma },
            { -M_PI, Dim::CYCLIC, M_PI, Dim::CYCLIC, num_angle_bins, dihedral_smoothing_sigma }
          }));
          break;
        case 3:
          features[std::make_pair(first_group, second_group)]
                = std::unique_ptr<cHistogram>(new cSmoothedHistogram({
            { 0, Dim::CLOSED, residue_interaction_cutoff, Dim::OPEN, num_dist_bins, displacement_smoothing_sigma },
            { 0, Dim::CLOSED, M_PI, Dim::CLOSED, num_angle_bins, theta_smoothing_sigma },
            { 0, Dim::CLOSED, M_PI, Dim::CLOSED, num_angle_bins, theta_smoothing_sigma }
          }));
          break;
        case 4:
          features[std::make_pair(first_group, second_group)]
                = std::unique_ptr<cHistogram>(new cSmoothedHistogram({
            { 0, Dim::CLOSED, residue_interaction_cutoff, Dim::OPEN, num_dist_bins, displacement_smoothing_sigma },
            { -M_PI, Dim::CYCLIC, M_PI, Dim::CYCLIC, num_angle_bins, dihedral_smoothing_sigma },
            { 0, Dim::CLOSED, M_PI, Dim::CLOSED, num_angle_bins, theta_smoothing_sigma },
            { 0, Dim::CLOSED, M_PI, Dim::CLOSED, num_angle_bins, theta_smoothing_sigma }
          }));
          break;
        case 5:
          features[std::make_pair(first_group, second_group)]
                = std::unique_ptr<cHistogram>(new cSmoothedHistogram({
            { 0, Dim::CLOSED, residue_interaction_cutoff, Dim::OPEN, num_dist_bins, displacement_smoothing_sigma },
            { 0, Dim::CLOSED, M_PI, Dim::CLOSED, num_angle_bins, theta_smoothing_sigma },
            { 0, Dim::CLOSED, M_PI, Dim::CLOSED, num_angle_bins, theta_smoothing_sigma },
            { 0, Dim::CYCLIC, 2 * M_PI, Dim::CYCLIC, num_angle_bins, phi_smoothing_sigma },
            { 0, Dim::CYCLIC, 2 * M_PI, Dim::CYCLIC, num_angle_bins, phi_smoothing_sigma }
          }));
          break;
        case 6:
          features[std::make_pair(first_group, second_group)]
                = std::unique_ptr<cHistogram>(new cSmoothedHistogram({
            { 0, Dim::CLOSED, residue_interaction_cutoff, Dim::OPEN, num_dist_bins, displacement_smoothing_sigma },
            { -M_PI, Dim::CYCLIC, M_PI, Dim::CYCLIC, num_angle_bins, dihedral_smoothing_sigma },
            { 0, Dim::CLOSED, M_PI, Dim::CLOSED, num_angle_bins, theta_smoothing_sigma },
            { 0, Dim::CLOSED, M_PI, Dim::CLOSED, num_angle_bins, theta_smoothing_sigma },
            { 0, Dim::CYCLIC, 2 * M_PI, Dim::CYCLIC, num_angle_bins, phi_smoothing_sigma },
            { 0, Dim::CYCLIC, 2 * M_PI, Dim::CYCLIC, num_angle_bins, phi_smoothing_sigma }
          }));
          break;
        default:
          assert(false);
      }
    }
  }
}

void cResiduePairwiseHist::accumulate_histograms(cProtein *protein) {
  auto neighbour_residues = NeighbourResidues(protein, sidechain_cutoff_,
                                              static_cast<double>(num_dist_bins_ + 1) / num_dist_bins_);
  for (const auto &residue_neighbours : neighbour_residues) {
    const cResidue *first = residue_neighbours.first;
    for (const cResidue *second : residue_neighbours.second) {
      if (get_group(*first) <= get_group(*second))
        add_event(*first, *second, &get_histogram(*first, *second));
    }
  }
}

cResiduePairwiseHist::ResidueGroup cResiduePairwiseHist::get_group(const cResidue &residue) const {
  switch (residue.type()) {
    case cResidue::AMINO: {
      auto *res = dynamic_cast<const cAminoResidue *>(&residue);
      return ResidueGroup(residue.type(), res->aminoType());
    }
    default: {
      return ResidueGroup(residue.type(), 0);
    }
  }
}

double get_dihedral(const cVector3 &A,
                    const cVector3 &B,
                    const cVector3 &C,
                    const cVector3 &D) {
  double dihedral = 0.0;
  math_functions::computeDihedralAngle(A, B, C, D, &dihedral);
  return dihedral / 180 * M_PI;
}

double get_angle(const cVector3 &A,
                 const cVector3 &B,
                 const cVector3 &C) {
  double angle = 0.0;
  math_functions::computeAngle(A, B, C, &angle);
  return angle / 180 * M_PI;
}

double get_phi(const cVector3 &N,
               const cVector3 &CA,
               const cVector3 &C,
               const cVector3 &CB,
               const cVector3 &second_CA) {
  double phi = 0.0;
  cVector3 projection_on_the_plane = second_CA - CA;
  cVector3 n = (N - CA) ^ (C - CA);
  n.orthogonalize(&projection_on_the_plane);

  math_functions::computeAngle(N, CA, projection_on_the_plane - CA, &phi);
  if ((((N - CA) ^ (projection_on_the_plane)) | (CB - CA)) < 0.0) {
    phi = 360.0 - phi;
  }
  return phi / 180 * M_PI;
}


void get_N_CA_C_CB_positions(const cAminoResidue &residue,
                             cVector3 *N, cVector3 *CA, cVector3 *C, cVector3 *CB) {
  const cAtom *N_atom = residue.getAtomByLabel(cAminoResidue::N),
              *CA_atom = residue.getAtomByLabel(cAminoResidue::CA),
              *C_atom = residue.getAtomByLabel(cAminoResidue::C),
              *CB_atom = (residue.aminoType() == cAminoResidue::GLY)
                          ? residue.getAtomByLabel(cAminoResidue::HA1)
                          : residue.getAtomByLabel(cAminoResidue::CB);
  assert(N_atom);
  assert(CA_atom);
  assert(C_atom);
  assert(CB_atom);

  *N = N_atom->getPosition(),
  *CA = CA_atom->getPosition(),
  *C = C_atom->getPosition(),
  *CB = CB_atom->getPosition();
}

void add_event(const cAminoResidue &first,
               const cAminoResidue &second,
               bool skip_unordered,
               cHistogram *histogram) {
  cVector3 first_N, first_CA, first_C, first_CB;
  cVector3 second_N, second_CA, second_C, second_CB;
  get_N_CA_C_CB_positions(first, &first_N, &first_CA, &first_C, &first_CB);
  get_N_CA_C_CB_positions(second, &second_N, &second_CA, &second_C, &second_CB);

  double displacement = (first_CA - second_CA).norm();

  if (histogram->num_dims() == 1) {
    if (!skip_unordered || (&first) < (&second))
      histogram->add_event({ displacement });
    return;
  }

  double dihedral = get_dihedral(first_CB, first_CA, second_CA, second_CB);

  if (histogram->num_dims() == 2) {
    if (!skip_unordered || (&first) < (&second))
      histogram->add_event({ displacement, dihedral });
    return;
  }

  double first_theta = get_angle(first_CB, first_CA, second_CA);
  double second_theta = get_angle(second_CB, second_CA, first_CA);

  if (histogram->num_dims() == 3) {
    if (!skip_unordered || std::make_tuple(first_theta, &first)
                            < std::make_tuple(second_theta, &second))
      histogram->add_event({ displacement,
                             first_theta, second_theta });
    return;
  }

  if (histogram->num_dims() == 4) {
    if (!skip_unordered || std::make_tuple(first_theta, &first)
                            < std::make_tuple(second_theta, &second))
      histogram->add_event({ displacement, dihedral,
                             first_theta, second_theta });
    return;
  }

  double first_phi = get_phi(first_N, first_CA, first_C, first_CB, second_CA);
  double second_phi = get_phi(second_N, second_CA, second_C, second_CB, first_CA);

  if (histogram->num_dims() == 5) {
    if (!skip_unordered || std::make_tuple(first_theta, first_phi, &first)
                            < std::make_tuple(second_theta, second_phi, &second))
      histogram->add_event({ displacement,
                             first_theta, second_theta,
                             first_phi, second_phi });
    return;
  }

  if (histogram->num_dims() == 6) {
    if (!skip_unordered || std::make_tuple(first_theta, first_phi, &first)
                            < std::make_tuple(second_theta, second_phi, &second))
      histogram->add_event({ displacement, dihedral,
                             first_theta, second_theta,
                             first_phi, second_phi });
    return;
  }

  assert(false);
}

void add_event(const cAminoResidue &first,
               const cResidue &second,
               bool skip_unordered,
               cHistogram *histogram) {
  try {
    add_event(first, dynamic_cast<const cAminoResidue &>(second), skip_unordered, histogram);
  } catch (...) {
    // TODO: amino + not amino
  }
}

void cResiduePairwiseHist::add_event(const cResidue &first,
                                     const cResidue &second,
                                     cHistogram *histogram) {
  assert(get_group(first) <= get_group(second));
  if (first.chainId == second.chainId
      && std::abs(static_cast<int>(first.seqNumber) - static_cast<int>(second.seqNumber))
         <= static_cast<int>(skipped_neighbourhood_)) {
    return;  // skip
  }

  bool skip_unordered = get_group(first) == get_group(second);
  try {
    ::add_event(dynamic_cast<const cAminoResidue &>(first), second, skip_unordered, histogram);
  } catch (...) {
    try {
      ::add_event(dynamic_cast<const cAminoResidue &>(second), first, skip_unordered, histogram);
    } catch (...) {
      // TODO: not amino + not amino
    }
  }
}

cAtomPairwiseHist::cAtomPairwiseHist(bool residue_type_dependent, bool atom_place_dependent,
                                     size_t num_bins,
                                     double cutoff, size_t skipped_neighbourhood,
                                     double smoothing_sigma)
    : residue_type_dependent_(residue_type_dependent),
      atom_place_dependent_(atom_place_dependent),
      max_interaction_dist_(cutoff * (num_bins + 1) / num_bins),
      skipped_neighbourhood_(skipped_neighbourhood) {
  std::vector<std::pair<int, int>> groups;

  auto residue_types = get_residue_types();
  auto atom_types = get_atom_types();

  for (int residue_type : residue_types) {
    for (int atom_type : atom_types) {
      groups.push_back({ residue_type, atom_type });
    }
  }

  for (const auto &first_group : groups) {
    for (const auto &second_group : groups) {
      if (first_group <= second_group) {
        features[std::make_pair(first_group, second_group)]
              = std::unique_ptr<cHistogram>(new cSmoothedHistogram({
                  { 0, cSmoothedHistogram::Dimension::CLOSED,
                    cutoff, cSmoothedHistogram::Dimension::OPEN,
                    num_bins, smoothing_sigma * std::sqrt(2) }
                }));
      }
    }
  }
}

std::vector<int> get_all_amino_types() {
  std::vector<int> amino_types;

  for (int amino_type = 0; amino_type <= cAminoResidue::AMINO_TYPE_END; ++amino_type) {
    amino_types.push_back(amino_type);
  }
  return amino_types;
}

std::vector<int> cAtomPairwiseHist::get_residue_types() const {
  if (!residue_type_dependent_) {
    return { cAminoResidue::AMINO_TYPE_END };
  }
  return get_all_amino_types();
}

std::vector<int> cAtomPairwiseHist::get_atom_types() const {
  if (!atom_place_dependent_) {
    return { 'C', 'N', 'O', 'H', 'S', 'X' };
  }

  std::vector<int> atom_types;
  for (int atom_type = 0; atom_type <= cAminoResidue::ATOM_LABEL_END; ++atom_type) {
    atom_types.push_back(atom_type);
  }
  return atom_types;
}

void cAtomPairwiseHist::accumulate_histograms(cProtein *protein) {
  auto neighbours = FindIsolatedNeighbours(protein, max_interaction_dist_, skipped_neighbourhood_);

  for (const auto &atom_neighbours : neighbours) {
    const auto *atom = atom_neighbours.first;
    for (const auto *neighbour : atom_neighbours.second) {
      assert(atom != neighbour);
      if (std::make_pair(get_group(*atom), atom)
          <= std::make_pair(get_group(*neighbour), neighbour)) {
        add_event(*atom, *neighbour, &get_histogram(*atom, *neighbour));
      }
    }
  }
}

void cAtomPairwiseHist::add_event(const cAtom &first,
                                  const cAtom &second,
                                  cHistogram *histogram) {
  double displacement = (first.getPosition() - second.getPosition()).norm();
  histogram->add_event({ displacement });
}

char GetElementGroup(const cAtom &atom) {
  char element_symbol = atom.elementSymbol[0];
  switch (element_symbol) {
    case 'C':
      return element_symbol;
    case 'N':
      return element_symbol;
    case 'O':
      return element_symbol;
    case 'H':
      return element_symbol;
    case 'S':
      return element_symbol;
    default:
      return 'X';
  }
}

std::pair<int, int> cAtomPairwiseHist::get_group(const cAtom &atom) const {
  assert(atom.residue);
  const cResidue &residue = *atom.residue;

  int residue_type = cAminoResidue::AMINO_TYPE_END;
  if (residue_type_dependent_ && residue.type() == cResidue::AMINO) {
    residue_type = dynamic_cast<const cAminoResidue *>(&residue)->aminoType();
  }

  int atom_type = GetElementGroup(atom);
  if (atom_place_dependent_) {
    atom_type = (residue.type() == cResidue::AMINO)
                ? dynamic_cast<const cAminoResidue *>(&residue)->get_atom_label(atom)
                : cAminoResidue::ATOM_LABEL_END;
  }

  return std::make_pair(residue_type, atom_type);
}


cHydrogenBondsHist::cHydrogenBondsHist(size_t num_dist_bins, size_t num_angle_bins,
                                       double max_bond_length, size_t skipped_neighbourhood,
                                       double smoothing_sigma)
    : max_interaction_dist_(max_bond_length * (num_dist_bins + 1) / num_dist_bins),
      skipped_neighbourhood_(skipped_neighbourhood) {
  typedef cSmoothedHistogram::Dimension Dim;

  features[std::make_pair(0, 0)]
      = std::unique_ptr<cHistogram>(new cSmoothedHistogram({
          { 0, Dim::CLOSED, max_bond_length, Dim::OPEN, num_dist_bins, smoothing_sigma * std::sqrt(2) },
          { 0, Dim::CLOSED, M_PI, Dim::CLOSED, num_angle_bins, smoothing_sigma },
          { 0, Dim::CLOSED, M_PI, Dim::CLOSED, num_angle_bins, smoothing_sigma }
        }));
}

void cHydrogenBondsHist::add_event(const cAminoResidue &O_residue,
                                   const cAminoResidue &H_residue,
                                   cHistogram *histogram) {
  const cAtom *C_atom = O_residue.getAtomByLabel(cAminoResidue::C),
              *O_atom = O_residue.getAtomByLabel(cAminoResidue::O),
              *H_atom = H_residue.getAtomByLabel(cAminoResidue::H),
              *N_atom = H_residue.getAtomByLabel(cAminoResidue::N);
  assert(C_atom && O_atom && H_atom && N_atom);

  if (O_atom->chainId == H_atom->chainId
      && std::abs(static_cast<int>(O_residue.seqNumber)
                  - static_cast<int>(H_residue.seqNumber)) == 1
      && HaveCovalentBond(*C_atom, *N_atom))
    return;

  cVector3 C = C_atom->getPosition(),
           O = O_atom->getPosition(),
           H = H_atom->getPosition(),
           N = N_atom->getPosition();

  double bond_length = (O - H).norm();
  double donor_angle = get_angle(C, O, H);
  double acceptor_angle = get_angle(O, H, N);

  histogram->add_event({ bond_length, donor_angle, acceptor_angle });
}

void cHydrogenBondsHist::accumulate_histograms(cProtein *protein) {
  std::vector<cAtom *> O_atoms;
  std::vector<cAtom *> H_atoms;

  ForEachAmino(protein, [&O_atoms, &H_atoms](cAminoResidue &residue) {
    auto *atom_O = const_cast<cAtom *>(residue.getAtomByLabel(cAminoResidue::O));
    if (atom_O) {
      O_atoms.push_back(atom_O);
      atom_O->neighbours.clear();
    }
    if (residue.aminoType() == cAminoResidue::PRO)
      return;

    auto *atom_H = const_cast<cAtom *>(residue.getAtomByLabel(cAminoResidue::H));
    if (atom_H) {
      H_atoms.push_back(atom_H);
    }
  });

  cGrid::initializeNeighbours(H_atoms, O_atoms, max_interaction_dist_);

  for (const auto *O_atom : O_atoms) {
    for (const auto *H_atom : O_atom->neighbours) {
      if (O_atom->chainId == H_atom->chainId
          && std::abs(static_cast<int>(O_atom->residue->seqNumber)
                      - static_cast<int>(H_atom->residue->seqNumber))
             <= static_cast<int>(skipped_neighbourhood_))
        continue;

      const cAminoResidue *amino_O = dynamic_cast<const cAminoResidue *>(O_atom->residue);
      const cAminoResidue *amino_H = dynamic_cast<const cAminoResidue *>(H_atom->residue);
      assert(amino_O);
      assert(amino_H);
      add_event(*amino_O, *amino_H, &get_histogram(*amino_O, *amino_H));
    }
  }
}


cSolvationShellHist::cSolvationShellHist(double sa_radius,
                                         double sa_smoothed_margin,
                                         size_t num_dist_bins,
                                         size_t num_angle_bins,
                                         double cutoff,
                                         double smoothing_sigma)
    : sa_radius_(sa_radius), sa_smoothed_margin_(sa_smoothed_margin),
      max_interaction_dist_(cutoff * (num_dist_bins + 1) / num_dist_bins) {
  auto residue_types = get_all_amino_types();
  int water_type = -1;

  for (int residue_type : residue_types) {
    typedef cSmoothedHistogram::Dimension Dim;

    features[std::make_pair(water_type, residue_type)]
          = std::unique_ptr<cHistogram>(new cSmoothedHistogram({
              { sa_radius, Dim::CLOSED, cutoff, Dim::OPEN, num_dist_bins, smoothing_sigma * std::sqrt(2) },
              { 0, Dim::CLOSED, M_PI, Dim::CLOSED, num_angle_bins, smoothing_sigma },
            }));
  }
}

int cSolvationShellHist::get_group(const cResidue &residue) const {
  switch (residue.type()) {
    case cResidue::WATER:
      return -1;
    case cResidue::AMINO:
      return dynamic_cast<const cAminoResidue *>(&residue)->aminoType();
    default:
      return cAminoResidue::AMINO_TYPE_END;
  }
}

void cSolvationShellHist::accumulate_histograms(cProtein *protein) {
  auto water_residues = GetHydrationShell(protein,
                                          sa_radius_, sa_smoothed_margin_,
                                          max_interaction_dist_);
  std::vector<cAtom *> water_atoms;
  for (const auto &water : water_residues) {
    assert(water->getAtom("O"));
    water_atoms.push_back(water->getAtom("O"));
  }

  std::vector<cAtom *> CA_atoms;
  ForEachAmino(protein, [&CA_atoms](cAminoResidue &residue) {
    auto *atom_CA = const_cast<cAtom *>(residue.getAtomByLabel(cAminoResidue::CA));
    assert(atom_CA);
    atom_CA->neighbours.clear();
    CA_atoms.push_back(atom_CA);
  });

  cGrid::initializeNeighbours(water_atoms, CA_atoms, max_interaction_dist_);

  for (const auto CA_atom : CA_atoms) {
    cResidue *residue = CA_atom->residue;
    for (const cAtom *water_atom : CA_atom->neighbours) {
      cResidue *water_residue = water_atom->residue;
      add_event(*water_residue, *residue, &get_histogram(*water_residue, *residue));
    }
  }
}

void cSolvationShellHist::add_event(const cResidue &water,
                                    const cResidue &residue,
                                    cHistogram *histogram) {
  auto amino_residue = dynamic_cast<const cAminoResidue *>(&residue);
  assert(amino_residue);
  const cAtom *O_atom = const_cast<cResidue &>(water).getAtom("O"),
              *CA_atom = amino_residue->getAtomByLabel(cAminoResidue::CA),
              *CB_atom = (amino_residue->aminoType() == cAminoResidue::GLY)
                          ? amino_residue->getAtomByLabel(cAminoResidue::HA1)
                          : amino_residue->getAtomByLabel(cAminoResidue::CB);
  assert(CA_atom);
  assert(CB_atom);
  assert(O_atom);

  cVector3 CA = CA_atom->getPosition(),
           CB = CB_atom->getPosition(),
           O = O_atom->getPosition();

  double displacement = (CA - O).norm();
  double theta = get_angle(CA, CB, O);

  histogram->add_event({ displacement, theta }, O_atom->occupancy);
}
