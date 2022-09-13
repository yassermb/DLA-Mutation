/*************************************************************************\

Sergei Grudinin, 2012
Mikhail Karasikov, 2016
All Rights Reserved.

\**************************************************************************/

#include "cProtein.hpp"

#include <stdio.h> // for file operations
//#include "bessel.h" // for bessel functions
//#include "besselfrac.h" // for bessel functions
#include <assert.h>
#include <string.h>
#include <cstdlib>

#include "cGrid.hpp"
#include "cAlgorithmTemplates.hpp"
#include "cRotamerLibrary.hpp"
#include "mathFunctions.hpp"
#include "cSparseMatrixOutput.hpp"
#include "energyModel.hpp"


using std::vector;
using std::set;


cProtein::~cProtein() {
  for (size_t j = 0; j < chains_.size(); ++j)
    delete chains_[j];

  if (rotamerLibrary)
    delete rotamerLibrary;
}

cVector3 cProtein::getCM(size_t j, double *mass_) const {
  double mass = 0;
  cVector3 pos(0, 0, 0);
  for (const auto &atom : chains_[j]->atoms()) {
    mass += atom.getMass();
    pos += atom.getPosition() * atom.getMass();
  }
  if (mass_)
    *mass_ = mass;
  return pos / mass;
}

cVector3 cProtein::getCM(double *mass_) const {
  double mass = 0;
  cVector3 pos(0, 0, 0);
  for (const auto &atom : atoms()) {
    mass += atom.getMass();
    pos += atom.getPosition() * atom.getMass();
  }
  if (mass_)
    *mass_ = mass;
  return pos / mass;
}

bool cProtein::writePDB(const std::string &pdbName, bool reindex) const {
  FILE *f = fopen(pdbName.c_str(), "w");
  if (!f) {
    perror("File writing error");
    return false;
  }
  if (reindex){
	  int n = 1;
	  for (const auto &atom : atoms()) {
		  cAtom at = atom; 
		  at.atomSerial = n;
		  at.toFile(f);
		  n++;
	  }
  }
  else {

	  for (const auto &atom : atoms())
		  atom.toFile(f);
  }

  fclose(f);
  return true;
}

size_t cProtein::numResidues() const {
  size_t num_residues = 0;
  for (size_t j = 0; j < chains_.size(); ++j)
    num_residues += chains_[j]->numResidues();
  return num_residues;
}

size_t cProtein::numAtoms() const {
  size_t num_atoms = 0;
  for (size_t j = 0; j < chains_.size(); ++j)
    num_atoms += chains_[j]->numAtoms();
  return num_atoms;
}

bool cProtein::isBackboneValid() const {
  for (const cProteinChain *chain : chains_) {
    if (chain && !chain->isBackboneValid())
      return false;
  }
  return true;
}

bool cProtein::initRotamersFromDunbrackBDLib(const std::string &libname) {
  assert(!rotamerLibrary);
  rotamerLibrary = new cPDBDunbrackBDLib();
  if (!rotamerLibrary->initialize(libname)) {
    std::cerr << "Loading Dunbrack library <" + libname + "> error\n";
    if (errno)
      perror(NULL);
    return false;
  }
	for (size_t j = 0; j < chains_.size(); j++) // loop through protein's chains
    chains_[j]->handleBackbone();

	for (size_t i = 0; i < vExpectedPhiPsi.size(); ++i) {
    auto resIt = residues().begin();
    while (resIt->chainId != vExpectedPhiPsi[i].chainId
            || resIt->seqNumber != vExpectedPhiPsi[i].resSeqNumber)
      ++resIt;
    if (std::abs(resIt->getPhi() - vExpectedPhiPsi[i].phi) > 0.01) {
      std::cerr << "ERROR! Chain: " << resIt->chainId << ", residue: " << resIt->seqNumber
                << ": Incorrect PHI value: " << resIt->getPhi()
                << " instead of " << vExpectedPhiPsi[i].phi << std::endl;
      exit(1);
    }
    if (std::abs(resIt->getPsi() - vExpectedPhiPsi[i].psi) > 0.01) {
      std::cerr << "ERROR! Chain: " << resIt->chainId << ", residue: " << resIt->seqNumber
                << ": Incorrect PSI value: " << resIt->getPsi()
                << " instead of " << vExpectedPhiPsi[i].psi << std::endl;
      exit(1);
    }
  }

  for (auto &residue : residues()) {
    cAminoResidue *amino_residue = dynamic_cast<cAminoResidue *>(&residue);
    if (!amino_residue)
      continue;

    amino_residue->predictSideChain(*rotamerLibrary);
    amino_residue->fitResidueGeometry(); // optionally
  }
  return true;
}

void cProtein::complete(bool withHydrogens) {
  for (cProteinChain *chain : chains_) {
    if (chain)
      chain->complete(withHydrogens);
  }
}

bool cProtein::setRotamers(const std::string &rotamers_filename) {
  vector<double> rotamers;
  if (!getVectorFromFile(rotamers_filename, &rotamers)) {
    fprintf(stderr, "cProtein::setRotamers error : Cant read vector from %s\n",
                                                      rotamers_filename.c_str());
    return false;
  }
  vector<cAminoResidue *> residues;
  energy::initializeNeighbourResidues(this, &residues);
  if (residues.size() != rotamers.size()) {
    fprintf(stderr, "cProtein::setRotamers error : wrong number of rotamers\n"
                    "Perhaps you are trying to download the file %s rotamers"
                    "corresponding to the other molecule.\n",
                                                      rotamers_filename.c_str());
    return false;
  }
  for (size_t resIdx = 0; resIdx < residues.size(); ++resIdx) {
    size_t rotamerIdx = static_cast<size_t>(rotamers[resIdx]);
    if (residues[resIdx]->numRotamers() > rotamerIdx) {
      residues[resIdx]->setRotamer(rotamerIdx);
      continue;
    }
    fprintf(stderr, "cProtein::setRotamers error : file rotamers %s "
                    "is not compatible with this molecule.\n",
                                                      rotamers_filename.c_str());
    return false;
  }
  return true;
}

cProtein::resItGen cProtein::residues() {
  return cProtein::resItGen(
    chains_,
    [](const cProteinChain &chain) { return chain.residues(); }
  );
}

const cProtein::resItGen cProtein::residues() const {
  return const_cast<cProtein &>(*this).residues();
}

cProtein::atomItGen cProtein::atoms() {
  return cProtein::atomItGen(
    chains_,
    [](const cProteinChain &chain) { return chain.atoms(); }
  );
}

const cProtein::atomItGen cProtein::atoms() const {
  return const_cast<cProtein &>(*this).atoms();
}
