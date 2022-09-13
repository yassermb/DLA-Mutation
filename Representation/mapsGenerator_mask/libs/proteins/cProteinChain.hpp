#pragma once

#include "cResidue.hpp"


class cProteinChain {
 public:
  explicit cProteinChain(char chainId_) : chainId(chainId_) {}

  cProteinChain(char chainId_, const std::vector<cResidue *> _residues)
    : chainId(chainId_), residues_(_residues) {}

  ~cProteinChain() {
    for (const auto &residue : residues_) {
      if (residue)
        delete residue;
    }
  }

  void addResidue(cResidue *residue);

  size_t numAtoms() const; // slow function O(numResidues)
  size_t numResidues() const { return residues_.size(); }

  void handleBackbone(); // Handle dihedral angles phi and psi

  void complete(bool withHydrogens = false);
  bool isBackboneValid() const;

  const char chainId;

  typedef FirstLevelIteratorGenerator<cResidue> resItGen;
  typedef HighLevelIteratorGenerator<cResidue, cAtom, cResidue::atomItGen> atomItGen;
  resItGen residues();
  const resItGen residues() const;
  atomItGen atoms();
  const atomItGen atoms() const;

 private:
  std::vector<cResidue *> residues_;
  //TODO: MODEL, Segment
};
