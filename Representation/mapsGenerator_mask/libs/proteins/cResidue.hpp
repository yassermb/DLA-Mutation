#pragma once

#include <map>

#include "cAtom.hpp"
#include "hierarchical_iterator.hpp"
#include "cSpatialTransform.hpp"

class cRotamerLibrary;
class cAminoResidue;
class cPDBNucleicResidue;
class cPDBIonResidue;
class cPDBWaterResidue;
class cPDBOtherResidue;

class cResidue {
 public:
  enum eType { AMINO = 0, NUCLEIC, WATER, ION, OTHER };

  cResidue(const std::string &resName, char chainId,
              const std::string &segId, size_t seqNumber, cResidue* previous);

  virtual ~cResidue() {
    for (size_t i = 0; i < atoms_.size(); ++i)
      delete atoms_[i];
  }

  cResidue(const cResidue &other) = delete;
  void operator=(const cResidue &other) = delete;

  static cResidue* create(const std::string &resName, char chainId,
                             const std::string &segId, size_t seqNumber, cResidue* previous);

  char getResCode() const;  // returns 'x' on failure
    size_t getResSeq() const {return seqNumber;}  

  cAtom* getAtom(const std::string &name) const ; // returns NULL on failure

  cAtom* getAtom(size_t i) const { return atoms_[i]; }

  virtual void addAtom(cAtom *atom);

  size_t numAtoms() const { return atoms_.size(); }

  virtual eType type() const = 0;

  virtual void handlePhiPsi(cResidue *nextRes) {}
  virtual void handlePhiPsi(cAminoResidue *nextRes) {}

  double getPhi() const { return phi; } // correct value in [-180, 180]
  double getPsi() const { return psi; } // 181.0 if value is not initialized
  const cResidue* getPrevious() const { return previous; }
  // iterators: atoms().begin(), atoms().end()
  typedef FirstLevelIteratorGenerator<cAtom> atomItGen;
  atomItGen atoms();
  const atomItGen atoms() const;

  size_t             residueIndex;

  const std::string  resName;    // ALA, ASN, HOH, etc...
  const char         chainId;    // Chain identifier
  const std::string  segId;      // Segment identifier, left-justified
  const size_t       seqNumber;  // Residue sequence number

    //function for oriented graphs
    cVector3 origin;
    cVector3 vX, vY, vZ;
    cSpatialTransform localToWorld, worldToLocal;
    virtual  bool setLocalFrame(bool orient = true) {return false;}

 protected:
  const cResidue*    previous;   //
  double phi; // 181.0 if not initialized
  double psi;

  virtual void removeAtom(std::vector<cAtom *>::iterator it);

  std::vector<cAtom *> atoms_;

 private:
  static const std::map<std::string, eType> mResType;
  static const std::map<std::string, char> mResCode;
};


class cOtherResidue : public cResidue {
 public:
  cOtherResidue(const std::string &resName, char chainId,
                const std::string &segId, size_t seqNumber, cResidue* previous)
    : cResidue(resName, chainId, segId, seqNumber, previous) {}

  eType type() const { return eType::OTHER; }
};

class cIonResidue : public cResidue {
 public:
  cIonResidue(const std::string &resName, char chainId,
              const std::string &segId, size_t seqNumber, cResidue* previous)
    : cResidue(resName, chainId, segId, seqNumber, previous) {}

  eType type() const { return eType::ION; }
};

class cWaterResidue : public cResidue {
 public:
  cWaterResidue(const std::string &resName, char chainId,
                const std::string &segId, size_t seqNumber, cResidue* previous)
    : cResidue(resName, chainId, segId, seqNumber, previous) {}

  eType type() const { return eType::WATER; }
};

class cNucleicResidue : public cResidue {
 public:
  cNucleicResidue(const std::string &resName, char chainId,
                  const std::string &segId, size_t seqNumber, cResidue* previous)
    : cResidue(resName, chainId, segId, seqNumber, previous) {}

  eType type() const { return eType::NUCLEIC; }
};
