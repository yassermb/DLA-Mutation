/*************************************************************************\

 Sergei Grudinin, NANO-D, 2010
 All Rights Reserved.

 \**************************************************************************/

#pragma once

#include <set>
#include <vector>
#include <sstream>

#include "cVector3.hpp"


class cResidue;

class cAtom {
 public:
  cAtom() : residue(NULL), charmmMass_(-1) {}
  cAtom(cAtom* atom);

  const cVector3& getPosition() const { return position; }
  void getPosition(double *x, double *y, double *z) const { position.get(x, y, z); }
  void setPosition(const cVector3 &p) { position=p; }
  void setPosition(double x, double y, double z) { position.set(x, y, z); }

  size_t        atomSerial;
  char          name[5]; // 4+'\0'
  char          altLoc;
  char          resName[4]; // 3+'\0'
  char          chainId;
  size_t        resSeqNumber;
  char          insertion;
  float         occupancy;
  float         bFactor;
  char          segId[5]; // 4+'\0'
  char          elementSymbol[3]; // 2+'\0'
  char          chargeOnAtom[3]; // 2+'\0'

  //int                                   rigidGroupNumber;

  int                                   getCharge() const;
  //float                                 getSolvationRadius() const { return solvationRadius; }
  //void                                  setSolvationRadius(float r) { solvationRadius = r; }
  float                                 getElementRadius() const { return elementRadius; }
  void                                  setElementRadius(float r) { elementRadius = r; }
  float                                 getMass() const { return mass; }
  void                                  setMass(float m) { mass = m; }

  friend std::ostream& operator<<(std::ostream & o, const cAtom &);

  void                                  toFile(FILE *f) const;
  bool                                  parseFromPDBField(const std::string &field);

  // some add info
  bool          isBackbone() const;
  bool          isSidechain() const;
  bool          isCA() const ;
  bool          isHydrogen() const;
  bool          isWaterOxygen() const;

  std::set<cAtom *> neighbours;

  cResidue *residue;

  // CHARMM
  double getCharmmMass() const { return charmmMass_; }
  double getCharmmMinEnergy() const { return charmmMinEnergy_; }
  double getCharmmMinRadius() const { return charmmMinRadius_; }
  double getCharmmPartCharge() const { return charmmPartCharge_; }

  void setCharmmMass(double charmmMass) { charmmMass_ = charmmMass; }
  void setCharmmMinEnergy(double charmmEnergy) { charmmMinEnergy_ = charmmEnergy; }
  void setCharmmMinRadius(double charmmRadius) { charmmMinRadius_ = charmmRadius; }
  void setCharmmPartCharge(double charmmPartCharge) { charmmPartCharge_ = charmmPartCharge; }

 private:
  cVector3 position;
  float elementRadius;
  //float solvationRadius;
  float mass;
  double charmmMass_;
  double charmmMinEnergy_;
  double charmmMinRadius_;
  double charmmPartCharge_;
};

// read string "ATOM  ...  \n"
// throws std::exception
std::istream& operator>>(std::istream & i, cAtom &);
