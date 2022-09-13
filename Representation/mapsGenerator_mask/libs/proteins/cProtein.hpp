/*************************************************************************\

Sergei Grudinin, 2012
All Rights Reserved.

\**************************************************************************/


#pragma once
#include <string>

#include "cProteinChain.hpp"


class cRotamerLibrary;
class cClassicalParserPDB;

class cProtein {
  friend cClassicalParserPDB;

 public:
  cProtein() : rotamerLibrary(NULL) {}
  ~cProtein();

  size_t numChains() const { return chains_.size(); }
  size_t numResidues() const;
  size_t numAtoms() const;

  cProteinChain* getChain(size_t j) { return chains_[j]; }

  cVector3 getCM(size_t j, double *mass = NULL) const;
  cVector3 getCM(double *mass = NULL) const;

  void complete(bool withHydrogens = false);
  bool writePDB(const std::string &, bool reindex = false) const;

  bool isBackboneValid() const;

  bool initRotamersFromDunbrackBDLib(const std::string &libname);
  bool setRotamers(const std::string &);

  typedef HighLevelIteratorGenerator<cProteinChain, cResidue, cProteinChain::resItGen> resItGen;
  typedef HighLevelIteratorGenerator<cProteinChain, cAtom, cProteinChain::atomItGen> atomItGen;
  resItGen residues();
  const resItGen residues() const;
  atomItGen atoms();
  const atomItGen atoms() const;

 private:
  std::vector<cProteinChain *> chains_;

  struct expectedPhiPsiContext {
    char chainId;
    size_t resSeqNumber;
    double phi;
    double psi;
  };
  std::vector<expectedPhiPsiContext> vExpectedPhiPsi;

  cRotamerLibrary *rotamerLibrary;
};
