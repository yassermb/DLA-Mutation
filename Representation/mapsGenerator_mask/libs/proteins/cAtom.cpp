/*************************************************************************\

 Sergei Grudinin, NANO-D, 2010
 All Rights Reserved.

 \**************************************************************************/

#include "cAtom.hpp"

#include <string.h>
#include <stdio.h>
#include <sstream>
#include <cmath>

#include "cPeriodicTable.hpp"
#include "cResidue.hpp"

const float kDefaultAtomRadius = 1.0;
const float kDefaultAtomMass = 1.0;


cAtom::cAtom(cAtom* atom) {
  atomSerial = atom->atomSerial;
  strcpy(name, atom->name);
  strcpy(resName, atom->resName);
  resSeqNumber = atom->resSeqNumber;
  chainId = atom->chainId;
  strcpy(segId, atom->segId);
  insertion = atom->insertion;
  altLoc = atom->altLoc;
  strcpy(elementSymbol, atom->elementSymbol);
  strcpy(chargeOnAtom, atom->chargeOnAtom);
  occupancy = atom->occupancy;
  bFactor = atom->bFactor;
  elementRadius = atom->elementRadius;
  setPosition (atom->getPosition());
}

bool cAtom::isCA() const {
  return !strcmp(name, "CA");
}

bool cAtom::isBackbone() const {
  int name_code = 0, i = 0;
  while (name[i])
    name_code = name_code << 8 | name[i++];

  switch (name_code) {
    case 'C' << 8 | 'A':
      return true;
    case 'C':
      return true;
    case 'N':
      return true;
    case 'O':
      return true;
    case ('O' << 8 | 'X') << 8 | 'T':
      return true;
    case 'H' << 8 | 'N':
      return true;
    case 'H':
      return true;
    case 'H' << 8 | 'A':
      return true;
    default:
      return false;
  }
  // if (!strcmp(name, "CA"))
    // return true;

  // if (!strcmp(name, "C"))
    // return true;

  // if (!strcmp(name, "N"))
    // return true;

  // if (!strcmp(name, "O"))
    // return true;

  // //if (!strncmp(name, "OT",2))
  // if (!strcmp(name, "OXT"))
    // return true;

  // if (!strcmp(name, "HN"))
    // return true;

  // if (!strcmp(name, "H"))
    // return true;

  // if (!strcmp(name, "HA"))
    // return true;

  // return false;
}

bool cAtom::isSidechain() const {
  return !isBackbone();
}

bool cAtom::isHydrogen() const {
  return !strcmp(elementSymbol, "H");
}

bool cAtom::isWaterOxygen() const {
  return residue->type() == cResidue::eType::WATER
         && !strcmp(elementSymbol, "O");
}

std::ostream& operator<<(std::ostream &out, const cAtom &atom) {
  double x,y,z;
  atom.getPosition(&x, &y, &z);
  char buffer[82];
  size_t len = strlen(atom.name);

  if (len == 4) {
    sprintf(buffer, "ATOM  %5zd %-4s%c%3s %c%4zd%c   %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s%2s%2s\n",
      atom.atomSerial, atom.name, atom.altLoc, atom.resName, atom.chainId,
      atom.resSeqNumber, atom.insertion, x, y, z, atom.occupancy, atom.bFactor,
      atom.segId, atom.elementSymbol, atom.chargeOnAtom);
  } else {
    sprintf(buffer, "ATOM  %5zd  %-3s%c%3s %c%4zd%c   %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s%2s%2s\n",
      atom.atomSerial, atom.name, atom.altLoc, atom.resName, atom.chainId,
      atom.resSeqNumber, atom.insertion, x, y, z, atom.occupancy, atom.bFactor,
      atom.segId, atom.elementSymbol, atom.chargeOnAtom);
  }
  out << buffer;
  return out;
};

bool cAtom::parseFromPDBField(const std::string &field_) {
  std::string field(field_);
  field.resize(81, ' ');
  field[80] = '\0';
  name[0] = resName[0] = segId[0] = elementSymbol[0] = chargeOnAtom[0] = '\0';

  sscanf(field.data() + 78, "%2s", chargeOnAtom);
  field[78] = '\0';
  if (sscanf(field.data() + 76, "%2s", elementSymbol) == EOF ||  cPeriodicTable::getMass(elementSymbol)==-1 ) {
    char element[] = { field.data()[12], field.data()[13], 0 };
    sscanf(element, "%2s", elementSymbol);
  }
  field[76] = '\0';
  sscanf(field.data() + 72, "%4s", segId);
  field[72] = '\0';
  sscanf(field.data() + 60, "%6f", &bFactor);
  field[60] = '\0';
  sscanf(field.data() + 54, "%6f", &occupancy);
  field[54] = '\0';
  double x, y, z;
  if (sscanf(field.data() + 46, "%8lf", &z) == EOF || std::isnan(z))
    return false;
  field[46] = '\0';
  if (sscanf(field.data() + 38, "%8lf", &y) == EOF || std::isnan(y))
    return false;
  field[38] = '\0';
  if (sscanf(field.data() + 30, "%8lf", &x) == EOF || std::isnan(x))
    return false;
  field[30] = '\0';
  setPosition(x, y, z);
  insertion = field[26];
  field[26] = '\0';
  if (sscanf(field.data() + 22, "%4zd", &resSeqNumber) == EOF)
    return false;
  chainId = field[21];
  field[21] = '\0';
  if (sscanf(field.data() + 17, "%3s", resName) == EOF)
    return false;
  altLoc = field[16];
  field[16] = '\0';
  if (sscanf(field.data() + 12, "%4s", name) == EOF)
    return false;
  field[12] = '\0';
  if (sscanf(field.data() + 6, "%5zd", &atomSerial) == EOF)
    return false;

  float mass = cPeriodicTable::getMass(elementSymbol);
  if (mass < 0)
    mass = kDefaultAtomMass;
  setMass(mass);
  float elementRadius = cPeriodicTable::getRadius(elementSymbol);
  if (elementRadius < 0)
    elementRadius = kDefaultAtomRadius;
  setElementRadius(elementRadius);

  return true;
}

class AtomParsingException : public std::exception {
 public:
  AtomParsingException(const std::string &field) : field_(field) {};
  ~AtomParsingException() throw() {}
  const char* what() const throw() {
    std::string what_string = "AtomParsingException with field:\n";
    return (what_string + field_).c_str();
  }
 private:
  std::string field_;
};

std::istream& operator>>(std::istream &in, cAtom &atom) {
  std::string field;
  std::getline(in, field);
  if (!atom.parseFromPDBField(field))
    throw AtomParsingException(field);
  return in;
}

int cAtom::getCharge() const {
  int charge = 0;
  sscanf(chargeOnAtom, "%d", &charge);
  return charge;
}

void cAtom::toFile(FILE *f) const {
  std::stringstream buffer;
  buffer << *this;
  fputs(buffer.str().c_str(), f);
};
