#include "cProteinChain.hpp"

#include <assert.h>
#include <algorithm>

#include "cAminoResidue.hpp"
#include "mathFunctions.hpp"


const double kMaxCNBondLength = 1.791;
const double kMinCNBondLength = 1.189;


void cProteinChain::addResidue(cResidue *residue) {
  residues_.push_back(residue);
}

size_t cProteinChain::numAtoms() const {
  size_t atoms_number = 0;
  for (size_t i = 0; i < residues_.size(); ++i) {
    if (residues_[i])
      atoms_number += residues_[i]->numAtoms();
  }
  return atoms_number;
}

void cProteinChain::handleBackbone() {
  for (size_t i = 0; i < residues_.size() - 1; ++i) {
    if (residues_[i] && residues_[i + 1])
      residues_[i]->handlePhiPsi(residues_[i + 1]);
  }
}

void cProteinChain::complete(bool withHydrogens) {
  for (auto &residue : residues()) {
    cAminoResidue *amino_residue = dynamic_cast<cAminoResidue *>(&residue);
    if (amino_residue)
      amino_residue->complete(withHydrogens);
  }

  // Reconstruct O and H atoms in backbone
  for (size_t i = 0; i + 1 < residues_.size(); ++i) {
    auto first_res = residues_[i];
    auto second_res = residues_[i + 1];

    if (!first_res || !second_res
        || std::abs(static_cast<int>(first_res->seqNumber)
                    - static_cast<int>(second_res->seqNumber)) != 1
        || first_res->type() != cResidue::AMINO
        || second_res->type() != cResidue::AMINO)
      continue;

    cAminoResidue *first = dynamic_cast<cAminoResidue *>(first_res);
    cAminoResidue *second = dynamic_cast<cAminoResidue *>(second_res);
    assert(first && second);

    if (!first->isBackboneValid() || !second->isBackboneValid())
      continue;

    assert(first->getAtomByLabel(cAminoResidue::CA));
    assert(first->getAtomByLabel(cAminoResidue::C));
    assert(first->getAtomByLabel(cAminoResidue::N));

    assert(second->getAtomByLabel(cAminoResidue::CA));
    assert(second->getAtomByLabel(cAminoResidue::C));
    assert(second->getAtomByLabel(cAminoResidue::N));

    cVector3 first_CA_pos = first->getAtomByLabel(cAminoResidue::CA)->getPosition(),
             first_C_pos = first->getAtomByLabel(cAminoResidue::C)->getPosition(),
             first_N_pos = first->getAtomByLabel(cAminoResidue::N)->getPosition();

    cVector3 second_CA_pos = second->getAtomByLabel(cAminoResidue::CA)->getPosition(),
             second_C_pos = second->getAtomByLabel(cAminoResidue::C)->getPosition(),
             second_N_pos = second->getAtomByLabel(cAminoResidue::N)->getPosition();

    if (std::min((first_C_pos - second_N_pos).norm(),
                 (second_C_pos - first_N_pos).norm()) > kMaxCNBondLength) {
      continue; //  not neighbours
    }

    // CA - N - C
    cVector3 CA_pos, N_pos, C_pos;

    if ((first_C_pos - second_CA_pos).norm() < (second_C_pos - first_CA_pos).norm()) {
      std::swap(first, second);
      C_pos = first_C_pos;
      if ((first_N_pos - second_CA_pos).norm() < (second_N_pos - first_CA_pos).norm())
        continue; // Such neighbours can't exist, skip this pair...
      N_pos = second_N_pos;
      CA_pos = second_CA_pos;
    } else {
      C_pos = second_C_pos;
      if ((first_N_pos - second_CA_pos).norm() > (second_N_pos - first_CA_pos).norm())
        continue; // Such neighbours can't exist, skip this pair...
      N_pos = first_N_pos;
      CA_pos = first_CA_pos;
    }

    if ((N_pos - C_pos).norm() < kMinCNBondLength || (N_pos - C_pos).norm() > kMaxCNBondLength)
      continue;

    if (!second->getAtomByLabel(cAminoResidue::O)) {
      cVector3 O_pos;
      math_functions::setDihedral(CA_pos, N_pos, C_pos, &O_pos, 1.233, 122.6, 0.0);

      second->addAtom(cAminoResidue::O, true);
      assert(second->getAtomByLabel(cAminoResidue::O));
      const_cast<cAtom*>(second->getAtomByLabel(cAminoResidue::O))->setPosition(O_pos);
    }

    if (withHydrogens
        && first->aminoType() != cAminoResidue::PRO // PRO doesn't have H atom
        && !first->getAtomByLabel(cAminoResidue::H)) {
      cVector3 H_pos;
      math_functions::setDihedral(CA_pos, C_pos, N_pos, &H_pos, 1.09, 120.0, 180.0);

      first->addAtom(cAminoResidue::H, true);
      assert(first->getAtomByLabel(cAminoResidue::H));
      const_cast<cAtom*>(first->getAtomByLabel(cAminoResidue::H))->setPosition(H_pos);
    }
  }
}

bool cProteinChain::isBackboneValid() const {
  for (const auto &residue : residues()) {
    auto *res = dynamic_cast<const cAminoResidue *>(&residue);
    if (res && !res->isBackboneValid())
      return false;
  }
  return true;
}

cProteinChain::resItGen cProteinChain::residues() {
  return cProteinChain::resItGen(residues_);
}

const cProteinChain::resItGen cProteinChain::residues() const {
  return const_cast<cProteinChain *>(this)->residues();
}

cProteinChain::atomItGen cProteinChain::atoms() {
  return cProteinChain::atomItGen(
    residues_,
    [](const cResidue &res) { return res.atoms(); }
  );
}

const cProteinChain::atomItGen cProteinChain::atoms() const {
  return const_cast<cProteinChain *>(this)->atoms();
}
