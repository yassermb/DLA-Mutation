#include "energyModel.hpp"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <algorithm>

#include "mathFunctions.hpp"
#include "cGrid.hpp"
#include "cAminoResidue.hpp"

const double kCutoffDistance = 5.0; // Angstroms

using std::set;
using std::vector;


namespace energy {

double vdwRadii(const char element[3]) {
  if (element[1]) {
    fprintf(stderr, "Amino acid can't contain element: %s\n", element);
    exit(1);
  }
  switch (element[0]) {
    case 'C':
      return 1.68;
    case 'N':
      return 1.28;
    case 'O':
      return 1.29;
    case 'H':
      return 0.49;
    case 'S':
      return 1.71;
    default:
      fprintf(stderr, "Unknown element: %s\n", element);
      exit(1);
  }
}

double vdwEnergy(const cAtom &atom_1, const cAtom &atom_2) {
  if (atom_1.getCharmmMass() < 0) {
    std::cerr << "CHARMM params was not initialized for atom\n" << atom_1
              << "vdw Energy was set to 0.0\n";
    return 0.0;
  }
  if (atom_2.getCharmmMass() < 0) {
    std::cerr << "CHARMM params was not initialized for atom\n" << atom_2
              << "vdw Energy was set to 0.0\n";
    return 0.0;
  }
  double d = (atom_1.getPosition() - atom_2.getPosition()).norm();
  double threshold = d / (vdwRadii(atom_1.elementSymbol)
                        + vdwRadii(atom_2.elementSymbol));
  if (threshold > 4 / 3)
    return 0.0;
  double E_ij = sqrt(atom_1.getCharmmMinEnergy() * atom_2.getCharmmMinEnergy());
  if (threshold > 10 / 9)
    return (E_ij / 4) * pow(9 * threshold - 10, 2) - E_ij;
  if (threshold > 1.0)
    return E_ij * pow(10 - 9 * threshold, 57.273 / (9 * E_ij)) - E_ij;
  if (threshold > 0.8254)
    return 57.273 * (1 - threshold);
  return 10.0;
}

double disulfideEnergy(const cAtom &atom_1, const cAtom &atom_2) {
  if (strcmp(atom_1.name, "SG") || strcmp(atom_2.name, "SG"))
    return 0.0;

  const cAminoResidue *cys_1 = dynamic_cast<cAminoResidue *>(atom_1.residue);
  const cAminoResidue *cys_2 = dynamic_cast<cAminoResidue *>(atom_2.residue);

  if (!cys_1 || !cys_2)
    return 0.0;

  const cAtom *SG_1 = cys_1->getAtomByLabel(cAminoResidue::eAtomLabel::SG),
                 *CB_1 = cys_1->getAtomByLabel(cAminoResidue::eAtomLabel::CB),
                 *SG_2 = cys_2->getAtomByLabel(cAminoResidue::eAtomLabel::SG),
                 *CB_2 = cys_2->getAtomByLabel(cAminoResidue::eAtomLabel::CB);
  double A_1, A_2, chi_3;
  bool ret = math_functions::computeAngle(SG_1->getPosition(),
                                          SG_2->getPosition(),
                                          CB_2->getPosition(), &A_1);
  assert(ret);
  ret = math_functions::computeAngle(SG_2->getPosition(),
                                     SG_1->getPosition(),
                                     CB_1->getPosition(), &A_2);
  assert(ret);
  ret = math_functions::computeDihedralAngle(CB_1->getPosition(),
                                             SG_1->getPosition(),
                                             SG_2->getPosition(),
                                             CB_2->getPosition(), &chi_3);
  assert(ret);
  double d = (SG_1->getPosition() - SG_2->getPosition()).norm();
  return std::min(0.0, 6 * (fabs(d - 2.06) + (fabs(A_1 - 105) + fabs(A_2 - 105)) / 100
                          + fabs(fabs(chi_3) - 90) / 140) - 11.4);
}

inline double cosAngle(double angle) {
  return cos(angle * M_PI / 180.0);
}

double hydrogenBondEnergy(const cAtom &atom_1, const cAtom &atom_2) {
  return 0.0;
/*
  if (atom_1.element[0] != 'O' || atom_2.element[0] != 'O' ||
      (atom_1.getPosition() - atom_2.getPosition()).norm() > 3.2)
    continue

  const cAminoResidue *hydroxyl
    = dynamic_cast<const cAminoResidue *>(atom_1.residue);
  const cAminoResidue *carboxyl
    = dynamic_cast<const cAminoResidue *>(atom_2.residue);

  const cAtom *C = carboxyl.getAtomByLabel(eAtomLabel::C);
  const cAtom *O = *pNeighbourIt;
  const cAtom *D = &(*atomIt);
  const cAtom *B;
  switch (hydroxyl.aminoType()) {
    case cAminoResidue::eAminoType::SER:
      B = hydroxyl.getAtomByLabel(eAtomLabel::CB);
      break;
    case cAminoResidue::eAminoType::THR:
      B = hydroxyl.getAtomByLabel(eAtomLabel::CB);
      break;
    case cAminoResidue::eAminoType::TYR:
      B = hydroxyl.getAtomByLabel(eAtomLabel::CZ);
      break;
    default:
      continue;
  }
// from RASP
// B and D -- hydroxyl, C and O -- carboxyl.
double hydrogenEnergy(const cAtom &B, const cAtom &D,
                      const cAtom &O, const cAtom &C) {
#ifdef HYDROGEN_ENERGY_ON
      if ((*pNeighbourIt)->element[0] != 'O' || atomIt->element[0] != 'O' ||
          (atomIt->getPosition() - (*pNeighbourIt)->getPosition()).norm() > 3.2)
        continue

      const cAminoResidue *hydroxyl
        = dynamic_cast<const cAminoResidue *>(&res);
      const cAminoResidue *carboxyl
        = dynamic_cast<const cAminoResidue *>((*pNeighbourIt)->residue);

      const cAtom *C = carboxyl.getAtomByLabel(eAtomLabel::C);
      const cAtom *O = *pNeighbourIt;
      const cAtom *D = &(*atomIt);
      const cAtom *B;
      switch (hydroxyl.aminoType()) {
        case cAminoResidue::eAminoType::SER:
          B = hydroxyl.getAtomByLabel(eAtomLabel::CB);
          break;
        case cAminoResidue::eAminoType::THR:
          B = hydroxyl.getAtomByLabel(eAtomLabel::CB);
          break;
        case cAminoResidue::eAminoType::TYR:
          B = hydroxyl.getAtomByLabel(eAtomLabel::CZ);
          break;
        default:
          continue;
      }
      energy += hydrogenEnergy(*B, *D, *O, *C);
#endif

  double alpha, beta;
  assert(math_functions::computeAngle(B.getPosition(),
                                      D.getPosition(),
                                      O.getPosition(), &alpha));
  assert(math_functions::computeAngle(C.getPosition(),
                                      O.getPosition(),
                                      D.getPosition(), &beta));
  return -1.8 * sqrt( (cosAngle(alpha - 111.5) - cosAngle(37.0))
                      * (cosAngle(beta - 120.0) - cosAngle(47.0))
                      / (1 - cosAngle(37.0))
                      / (1 - cosAngle(47.0)) );
*/
}

// side-chain --- frame
double frameEnergy(const cResidue &res) {
  double energy = 0.0;
  for (const auto &atom : res.atoms()) {
    if (atom.isBackbone())
      continue;

    for (const cAtom *neighbour : atom.neighbours) {
      energy += vdwEnergy(atom, *neighbour);
      energy += disulfideEnergy(atom, *neighbour);
      energy += hydrogenBondEnergy(atom, *neighbour);
    }
  }
  return energy;
}

// side-chain --- side-chain
double pairEnergy(const cResidue &res_1, const cResidue &res_2) {
  double energy = 0.0;
  for (const auto &first_atom : res_1.atoms()) {
    if (first_atom.isBackbone())
      continue;

    for (const auto &second_atom : res_2.atoms()) {
      if (second_atom.isBackbone() || &first_atom == &second_atom)
        continue;

      energy += vdwEnergy(first_atom, second_atom);
      energy += disulfideEnergy(first_atom, second_atom);
      energy += hydrogenBondEnergy(first_atom, second_atom);
    }
  }
  return energy;
}

double totalEnergy(const cResidue &res) {
  double energy = 0.0;
  for (const auto &atom : res.atoms()) {
    // if (atom.isBackbone())
      // continue;

    for (const cAtom *neighbour : atom.neighbours) {
      if (atom.isBackbone() && neighbour->isBackbone())
        continue;

      energy += vdwEnergy(atom, *neighbour);
      energy += disulfideEnergy(atom, *neighbour);
      energy += hydrogenBondEnergy(atom, *neighbour);
    }
  }
  return energy;
}

// RASP radii of the sidechain hemispheres
double scHemisphereRadii(cAminoResidue::eAminoType aminoType) {
  switch (aminoType) {
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
    default:
      fprintf(stderr, "Unknown amino-residue\n");
      exit(1);
  }
}

vector<vector<cAminoResidue *>>
    initializeNeighbourResidues(cProtein *protein, vector<cAminoResidue *> *residues) {
  residues->clear();
  // each residue is associated with CB atom: residue := CB->residue
  // so we keep only CB atoms instead of residues
  vector<cAtom *> CBAtoms;
  vector<cAtom *> backbone, sidechain;
  for (auto &atom : protein->atoms()) {
    atom.neighbours.clear();
    cAminoResidue *residue = dynamic_cast<cAminoResidue *>(atom.residue);
    // FYI: We calculate energy only with 20 standard amino residues
    if (!residue)
      continue;

    if (residue->numRotamers() < 2) {
      backbone.push_back(&atom);
      continue;
    }

    if (atom.isBackbone()) {
      backbone.push_back(&atom);
    } else {
      sidechain.push_back(&atom);
    }

    if (!strcmp(atom.name, "CB")) {
      CBAtoms.push_back(&atom);
      residues->push_back(residue);
    }
  }

  cGrid::initializeNeighbours(CBAtoms, CBAtoms,
    kCutoffDistance + 2 * scHemisphereRadii(cAminoResidue::eAminoType::ARG));

  vector<vector<cAminoResidue *>> neighbours(CBAtoms.size());
  for (size_t i = 0; i < CBAtoms.size(); ++i) {
    cAminoResidue *residue = residues->at(i);
    // RASP method of detecting residues that have effective contact
    for (set<cAtom *>::iterator it = CBAtoms[i]->neighbours.begin();
                                  it != CBAtoms[i]->neighbours.end(); ++it) {
      cAminoResidue *neighbour = dynamic_cast<cAminoResidue *>((*it)->residue);
      cVector3 CB = CBAtoms[i]->getPosition(),
               NbCB = (*it)->getPosition();
      if ((CB - NbCB).norm() > scHemisphereRadii(residue->aminoType())
                                + scHemisphereRadii(neighbour->aminoType())
                                + kCutoffDistance)
        continue;
#if 1
      neighbours[i].push_back(neighbour);
#else
      /* С Этим работает плохо. Необходимо перепроверить. */
      cVector3 CA = residue->getAtomByLabel(cAminoResidue::eAtomLabel::CA)->getPosition(),
               NbCA = neighbour->getAtomByLabel(cAminoResidue::eAtomLabel::CA)->getPosition();
      double alpha, beta;
      assert(math_functions::computeAngle(CA, CB, NbCB, &alpha));
      assert(math_functions::computeAngle(NbCA, NbCB, CB, &beta));
      if (alpha > 90 || beta > 90)
        neighbours[i].push_back(neighbour);
#endif
    }
    CBAtoms[i]->neighbours.clear();
  }
  cGrid::initializeNeighbours(backbone, sidechain, kCutoffDistance
                              + scHemisphereRadii(cAminoResidue::eAminoType::ARG));
  return neighbours;
}

void initializeNeighbourAtoms(cProtein *protein) {
  vector<cAtom *> atoms;
  for (auto &atom : protein->atoms()) {
    atom.neighbours.clear();
    atoms.push_back(&atom);
  }
  cGrid::initializeNeighbours(atoms, atoms, kCutoffDistance
                              + scHemisphereRadii(cAminoResidue::eAminoType::ARG)); // TUNE IT
}

size_t numClashes(const cProtein &protein) {
  vector<cAtom *> pAtoms;
  for (const auto &atom : protein.atoms()) {
    pAtoms.push_back(new cAtom(atom));
    pAtoms.back()->neighbours.clear();
  }

  cGrid::initializeNeighbours(pAtoms, pAtoms, kCutoffDistance);

  size_t clashes = 0;
  for (size_t i = 0; i < pAtoms.size(); ++i) {
    for (set<cAtom *>::iterator it = pAtoms[i]->neighbours.begin();
                                  it != pAtoms[i]->neighbours.end(); ++it) {
      if ((*it)->residue > pAtoms[i]->residue
          && (pAtoms[i]->isSidechain() || (*it)->isSidechain())
          && ((*it)->getPosition() - pAtoms[i]->getPosition()).norm()
              < 0.6 * (vdwRadii((*it)->elementSymbol)
                      + vdwRadii(pAtoms[i]->elementSymbol)))
        clashes++;
    }
  }
  for (size_t i = 0; i < pAtoms.size(); ++i)
    delete pAtoms[i];

  return clashes;
}


} // energy
