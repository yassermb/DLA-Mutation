#pragma once

#include "cResidue.hpp"


struct cRotamer {
  double chi[4]; // rotamer angles
  double prob; // rotamer probability
};

class cAminoResidue : public cResidue {
 public:
  enum eAminoType { ALA = 0, ARG, ASN, ASP, CYS, GLN, GLU, GLY, HIS, ILE,
                    LEU, LYS, MET, PHE, PRO, SER, THR, TRP, TYR, VAL, SEC, MSE,
                    AMINO_TYPE_END };

  cAminoResidue(const std::string &resName, char chainId,
                   const std::string &segId, size_t seqNumber, cResidue* previous_);

    virtual bool setLocalFrame(bool orient = true);

  eType type() const { return eType::AMINO; }

  eAminoType aminoType() const { return aminoType_; }

  void addAtom(cAtom *atom);
  // add all missing atoms in default state
  void complete(bool withHydrogens = false);

  void handlePhiPsi(cResidue *nextRes) { nextRes->handlePhiPsi(this); }
  void handlePhiPsi(cAminoResidue *nextRes);

  bool isBackboneValid() const;

  // initialize possible rotamers
  void predictSideChain(const cRotamerLibrary &);
  size_t numRotamers() const;
  // fit geometry of the residue from current rotamer
  void fitResidueGeometry();
  // |libEnergy| expresses the relative probability of the j-th
  // rotamer to the highest probability 0-th rotamer
  // return value: probability of the rotamer
  double setRotamer(size_t j, double *libEnergy = NULL);

  // get chi_1, chi_2, chi_3, or chi_4 dihedral angle
  double getChi(size_t n) const;

  static const std::map<std::string, eAminoType> mAminoType;

  enum eAtomLabel {
    C = 0, N, O, H, HT1, HT2, HT3,

    CA, HA, HA1, HA2, HA3,

    CB, HB, HB1, HB2, HB3,

    CG, CG1, CG2, HG, HG1, HG2, HG3,
    HG11, HG12, HG13,
    HG21, HG22, HG23,
    OG, OG1, SG,

    CD, CD1, CD2, HD, HD1, HD2, HD3,
    HD11, HD12, HD13, HD21, HD22, HD23,
    SD, OD1, OD2, ND1, ND2,

    CE, CE1, CE2, CE3, HE, HE1, HE2, HE3,
    HE21, HE22,
    OE1, OE2, NE, NE1, NE2,

    CZ, CZ2, CZ3, HZ, HZ1, HZ2, HZ3, NZ,

    CH2, CH3, NH1, NH2, OH, HH, HH11,
    HH12, HH2, HH21, HH22,
    HH31, HH32, HH33, O1, O2, OXT, SE,

    ATOM_LABEL_END
  };

  eAtomLabel get_atom_label(const cAtom &atom) const;

  // add missing atom in default state
  void addAtom(eAtomLabel atomLabel, bool withHydrogens = false);
  const cAtom* getAtomByLabel(eAtomLabel atomLabel) const;

 private:
  void removeAtom(std::vector<cAtom *>::iterator it);

  static const std::map<std::string, eAtomLabel> mAtomLabel;
  std::map<eAtomLabel, cAtom *> atom_map;

  eAminoType aminoType_;

  struct residueRebuildingContext {
    cAminoResidue::eAtomLabel A[4];
    double dist, angle, dihedral;
    size_t chi; // 1,2,3,4 or 0, if angle is known
  };

  struct cRotamersContext {
    ~cRotamersContext() { delete[] residueRebuildingData; }

    residueRebuildingContext *residueRebuildingData;
    size_t residueRebuildingDataSize;

    std::vector<cRotamer> rotamers;
  };

  cRotamersContext rotamersContext;
  std::map<eAtomLabel, eAtomLabel> mAtomAliases;

  static const residueRebuildingContext ALAresidueRebuildingData[];
  static const residueRebuildingContext ARGresidueRebuildingData[];
  static const residueRebuildingContext ASNresidueRebuildingData[];
  static const residueRebuildingContext ASPresidueRebuildingData[];
  static const residueRebuildingContext CYSresidueRebuildingData[];
  static const residueRebuildingContext GLNresidueRebuildingData[];
  static const residueRebuildingContext GLUresidueRebuildingData[];
  static const residueRebuildingContext GLYresidueRebuildingData[];
  static const residueRebuildingContext HISresidueRebuildingData[];
  static const residueRebuildingContext ILEresidueRebuildingData[];
  static const residueRebuildingContext LEUresidueRebuildingData[];
  static const residueRebuildingContext LYSresidueRebuildingData[];
  static const residueRebuildingContext METresidueRebuildingData[];
  static const residueRebuildingContext PHEresidueRebuildingData[];
  static const residueRebuildingContext PROresidueRebuildingData[];
  static const residueRebuildingContext SERresidueRebuildingData[];
  static const residueRebuildingContext THRresidueRebuildingData[];
  static const residueRebuildingContext TRPresidueRebuildingData[];
  static const residueRebuildingContext TYRresidueRebuildingData[];
  static const residueRebuildingContext VALresidueRebuildingData[];
  static const residueRebuildingContext SECresidueRebuildingData[];
  static const residueRebuildingContext MSEresidueRebuildingData[];
};
