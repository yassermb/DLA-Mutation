#include "NeighbourSearch.hpp"

#include <assert.h>
#include <algorithm>
#include <iterator>

#include "cGrid.hpp"


const double kMaxBondLength = 2.5;
const double kMinDistanceBetweenAtoms = 0.01;


double scHemisphereRadii(cAminoResidue::eAminoType aminoType) {
  switch (aminoType) {
    case cAminoResidue::eAminoType::ALA:
      // TODO: measure
      return 2.8;
    case cAminoResidue::eAminoType::ARG:
      return 6.3;
    case cAminoResidue::eAminoType::ASN:
      return 3.8;
    case cAminoResidue::eAminoType::ASP:
      return 3.8;
    case cAminoResidue::eAminoType::CYS:
      return 3.8;
    case cAminoResidue::eAminoType::GLN:
      return 4.1;
    case cAminoResidue::eAminoType::GLU:
      return 4.1;
    case cAminoResidue::eAminoType::GLY:
    // TODO: measure
      return 0.0;
    case cAminoResidue::eAminoType::HIS:
      return 3.8;
    case cAminoResidue::eAminoType::ILE:
      return 3.8;
    case cAminoResidue::eAminoType::LEU:
      return 3.8;
    case cAminoResidue::eAminoType::LYS:
      return 5.4;
    case cAminoResidue::eAminoType::MET:
      return 4.5;
    case cAminoResidue::eAminoType::PHE:
      return 4.5;
    case cAminoResidue::eAminoType::PRO:
      return 3.8;
    case cAminoResidue::eAminoType::SER:
      return 3.8;
    case cAminoResidue::eAminoType::THR:
      return 3.8;
    case cAminoResidue::eAminoType::TRP:
      return 5.6;
    case cAminoResidue::eAminoType::TYR:
      return 5.9;
    case cAminoResidue::eAminoType::VAL:
      return 3.8;
    case cAminoResidue::eAminoType::SEC:
      return 3.8;
    case cAminoResidue::eAminoType::MSE:
      return 4.5;
    default:
      assert(false && "Unknown amino-residue");
  }
  exit(1);
}

double contact_dist(cAminoResidue::eAminoType first,
                    cAminoResidue::eAminoType second) {
  return scHemisphereRadii(first) + scHemisphereRadii(second);
}

std::map<cResidue *, std::vector<cResidue *>> NeighbourResidues(cProtein *protein,
                                                                double sidechain_cutoff,
                                                                double residue_dist_multiplier) {
  assert(residue_dist_multiplier >= 1);

  std::map<cResidue *, std::vector<cResidue *>> neighbours;

  std::vector<cAtom *> representatives;

  ForEachAmino(protein, [&representatives](cAminoResidue &residue) {
    auto *atom_CA = const_cast<cAtom *>(residue.getAtomByLabel(cAminoResidue::CA));
    assert(atom_CA);
    atom_CA->neighbours.clear();
    representatives.push_back(atom_CA);
  });

  cGrid::initializeNeighbours(representatives, representatives,
    (sidechain_cutoff + contact_dist(cAminoResidue::ARG,
                                     cAminoResidue::ARG)) * residue_dist_multiplier);

  for (const auto representative : representatives) {
    cResidue *residue = representative->residue;
    // RASP method of detecting residues that have effective contact
    for (const cAtom *atom : representative->neighbours) {
      cResidue *neighbour = atom->residue;
      cVector3 first = representative->getPosition(),
               second = atom->getPosition();
      if ((sidechain_cutoff + contact_dist(dynamic_cast<cAminoResidue *>(residue)->aminoType(),
                                           dynamic_cast<cAminoResidue *>(neighbour)->aminoType()))
           * residue_dist_multiplier > (first - second).norm()) {
        neighbours[residue].push_back(neighbour);
      }
    }
    representative->neighbours.clear();
  }
  return neighbours;
}


// checks if there exists a covalent bond between two atoms,
// based on the distance between atoms, their radii and their specific properties
// co-developed with S. Artemova
bool HaveCovalentBond(const cAtom &first, const cAtom &second) {
  double distance2 = (first.getPosition() - second.getPosition()).norm2();

  if (first.residue->type() == cResidue::WATER
      && second.residue->type() == cResidue::WATER)
    return distance2 < 1.21;

  // then, we don't connect water with anything else
  if ((first.residue->type() == cResidue::WATER)
      ^ (second.residue->type() == cResidue::WATER))
    return false;

  // another rule, we don't connect things from different chains
  if (first.chainId != second.chainId)
    return false;

  if (first.isHydrogen() && second.isHydrogen())
    return false;

  if ((first.isHydrogen() || second.isHydrogen())
      && distance2 >= 1.21) {
    return false;
  }

  // we don't create S-S bridges
  if (!std::strncmp(first.name, "SG", 2) && !std::strncmp(second.name, "SG", 2))
    return false;

  // we don't connect things that are far in the sequence
  if (std::abs(static_cast<int>(first.resSeqNumber)
               - static_cast<int>(second.resSeqNumber)) > 1)
    return false;

  // here we don't want to have a connection between two "H" in a water
  // but there are some forcefields where it can be possible

  double distanceMax2 = (first.getElementRadius() + second.getElementRadius()) * 0.6;
  distanceMax2 *= distanceMax2;

  // we don't consider pairs which are too close also because that could be a mistake
  return distance2 <= distanceMax2 && distance2 > kMinDistanceBetweenAtoms;
}


// This function is developed to compute atoms in the topological distance |order|.
std::map<cAtom *, std::set<cAtom *>> FindTopologicalNeighbours(cProtein *protein,
                                                               size_t order) {
  std::vector<cAtom *> atoms;

  for (cAtom &atom : protein->atoms()) {
    atom.neighbours.clear();
    atoms.push_back(&atom);
  }

  cGrid::initializeNeighbours(atoms, atoms, kMaxBondLength);

  for (cAtom *atom : atoms) {
    std::set<cAtom *> not_covalent_neighbours;
    for (cAtom *neighbour : atom->neighbours) {
      if (!HaveCovalentBond(*atom, *neighbour))
        not_covalent_neighbours.insert(neighbour);
    }
    atom->neighbours = SetDifference(atom->neighbours, not_covalent_neighbours);
  }

  std::vector<std::vector<std::set<cAtom *>>> neighbours(
    order + 1,
    std::vector<std::set<cAtom *>>(atoms.size())
  );

  for (size_t i = 0; i < atoms.size(); ++i) {
    neighbours[0][i].insert(atoms[i]);
  }

  for (size_t current_order = 0; current_order < order; ++current_order) {
    for (size_t i = 0; i < atoms.size(); ++i) {
      for (cAtom *from : neighbours[current_order][i]) {
        for (cAtom *to : from->neighbours) {
          neighbours[current_order + 1][i].insert(to);
        }
      }
    }
  }

  for (cAtom *atom : atoms) {
    atom->neighbours.clear();
  }

  std::map<cAtom *, std::set<cAtom *>> result;

  for (size_t i = 0; i < atoms.size(); ++i) {
    for (size_t current_order = 0; current_order <= order; ++current_order) {
      auto &reachable_atoms = result[atoms[i]];

      for (cAtom *atom : neighbours[current_order][i]) {
        reachable_atoms.insert(atom);
      }
    }
  }

  return result;
}


std::map<cAtom *, std::set<cAtom *>> FindIsolatedNeighbours(cProtein *protein,
                                                            double cutoff,
                                                            size_t order) {
  auto topological_neighbours = FindTopologicalNeighbours(protein, order);

  std::vector<cAtom *> atoms;
  for (cAtom &atom : protein->atoms()) {
    atom.neighbours.clear();
    atoms.push_back(&atom);
  }
  cGrid::initializeNeighbours(atoms, atoms, cutoff);

  std::map<cAtom *, std::set<cAtom *>> result;

  for (cAtom *atom : atoms) {
    result[atom] = SetDifference(atom->neighbours, topological_neighbours[atom]);
  }

  return result;
}
