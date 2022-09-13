#include "cRotamerLibrary.hpp"

#include <assert.h>
#include <cstring>
#include <algorithm>


cResidue::cResidue(const std::string &resName_, char chainId_,
                         const std::string &segId_, size_t seqNumber_, cResidue* previous_)
    : resName(resName_), chainId(chainId_),
      segId(segId_), seqNumber(seqNumber_), phi(181.0), psi(181.0), previous(previous_){}

cResidue* cResidue::create(const std::string &resName_, char chainId_,
                           const std::string &segId_, size_t seqNumber_, cResidue* previous_) {
  std::map<std::string, cResidue::eType>::const_iterator it = mResType.find(resName_);
  if (it != mResType.end()) {
    if (it->second == eType::AMINO)
      return new cAminoResidue(resName_, chainId_, segId_, seqNumber_, previous_);
    if (it->second == eType::NUCLEIC)
      return new cNucleicResidue(resName_, chainId_, segId_, seqNumber_, previous_);
    if (it->second == eType::WATER)
      return new cWaterResidue(resName_, chainId_, segId_, seqNumber_, previous_);
    if (it->second == eType::ION)
      return new cIonResidue(resName_, chainId_, segId_, seqNumber_, previous_);
  }
  return new cOtherResidue(resName_, chainId_, segId_, seqNumber_, previous_);
}

void cResidue::addAtom(cAtom *new_atom) {
  assert(new_atom);
  new_atom->residue = this;

  if (new_atom->altLoc == ' ') {
    atoms_.push_back(new_atom);
    return;
  }

  auto it = std::find_if(atoms_.begin(), atoms_.end(),
                          [new_atom](cAtom *existing) {
                            return !std::strcmp(existing->name, new_atom->name);
                          });
  if (it == atoms_.end()) {
    atoms_.push_back(new_atom);
    return;
  }

  if (new_atom->occupancy > (*it)->occupancy) {
    removeAtom(it);
    atoms_.push_back(new_atom);
  } else {
    delete new_atom;
  }
}

void cResidue::removeAtom(std::vector<cAtom *>::iterator it) {
  delete *it;
  atoms_.erase(it);
}

char cResidue::getResCode() const {
  std::map<std::string, char>::const_iterator it = mResCode.find(resName);
  if (it != mResCode.end())
    return it->second;
  return 'x';
}

cAtom* cResidue::getAtom(const std::string &name) const {
  for (size_t i = 0; i < atoms_.size(); ++i) {
    if (name == atoms_[i]->name)
      return atoms_[i];
  }
  return NULL;
}

cResidue::atomItGen cResidue::atoms() {
  return cResidue::atomItGen(atoms_);
}

const cResidue::atomItGen cResidue::atoms() const {
  return cResidue::atomItGen(atoms_);
}

const std::map<std::string, cResidue::eType> cResidue::mResType = {

  {"CA", cResidue::eType::ION},
  {"NA", cResidue::eType::ION},
  {"CL", cResidue::eType::ION},
  {"K", cResidue::eType::ION},
  {"MG", cResidue::eType::ION},
  {"ZN", cResidue::eType::ION},

  {"HOH", cResidue::eType::WATER},
  {"HHO", cResidue::eType::WATER},
  {"OHH", cResidue::eType::WATER},
  {"H2O", cResidue::eType::WATER},
  {"OH2", cResidue::eType::WATER},
  {"WAT", cResidue::eType::WATER},
  {"TIP", cResidue::eType::WATER},
  {"TIP3", cResidue::eType::WATER},
  {"TIP4", cResidue::eType::WATER},
  {"TIP3P", cResidue::eType::WATER},
  {"TIP4P", cResidue::eType::WATER},
  {"SOL", cResidue::eType::WATER},

  {"A", cResidue::eType::NUCLEIC},
  {"+A", cResidue::eType::NUCLEIC},
  {"C", cResidue::eType::NUCLEIC},
  {"+C", cResidue::eType::NUCLEIC},
  {"G", cResidue::eType::NUCLEIC},
  {"+G", cResidue::eType::NUCLEIC},
  {"I", cResidue::eType::NUCLEIC},
  {"+I", cResidue::eType::NUCLEIC},
  {"T", cResidue::eType::NUCLEIC},
  {"+T", cResidue::eType::NUCLEIC},
  {"U", cResidue::eType::NUCLEIC},
  {"+U", cResidue::eType::NUCLEIC},

  {"ALA", cResidue::eType::AMINO},
  {"ARG", cResidue::eType::AMINO},
  {"ASN", cResidue::eType::AMINO},
  {"ASP", cResidue::eType::AMINO},
  {"CYS", cResidue::eType::AMINO},
  {"GLN", cResidue::eType::AMINO},
  {"GLU", cResidue::eType::AMINO},
  {"GLY", cResidue::eType::AMINO},
  {"HIS", cResidue::eType::AMINO},
  {"ILE", cResidue::eType::AMINO},
  {"LEU", cResidue::eType::AMINO},
  {"LYS", cResidue::eType::AMINO},
  {"MET", cResidue::eType::AMINO},
  {"PHE", cResidue::eType::AMINO},
  {"PRO", cResidue::eType::AMINO},
  {"SER", cResidue::eType::AMINO},
  {"THR", cResidue::eType::AMINO},
  {"TRP", cResidue::eType::AMINO},
  {"TYR", cResidue::eType::AMINO},
  {"VAL", cResidue::eType::AMINO},
  {"SEC", cResidue::eType::AMINO},
  {"MSE", cResidue::eType::AMINO}
};

const std::map<std::string, char> cResidue::mResCode = {

  {"ALA", 'A'},
  {"ARG", 'R'},
  {"ASN", 'N'},
  {"ASP", 'D'},
  {"CYS", 'C'},
  {"GLN", 'Q'},
  {"GLU", 'E'},
  {"GLY", 'G'},
  {"HIS", 'H'},
  {"ILE", 'I'},
  {"LEU", 'L'},
  {"LYS", 'K'},
  {"MET", 'M'},
  {"PHE", 'F'},
  {"PRO", 'P'},
  {"SER", 'S'},
  {"THR", 'T'},
  {"TRP", 'W'},
  {"TYR", 'Y'},
  {"VAL", 'V'}

};
