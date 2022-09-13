#pragma once

#include "cAminoResidue.hpp"


class cCharmmParams {
 public:
  enum eAtomName {
    H = 0,HC,HA,HT,LP,CT,C,CH1E,CH2E,CH3E,CR1E,CM,N,NR,NP,NH1E,NH2E,NH3E,
    NC2E,NH1,NH2,NH3,NC2,O,OC,OH1E,OH2E,OH1,OH2,OM,OT,OS,S,SH1E,CAL,FE,SE,DUMMY
  };

  struct toph19 {
    eAtomName charmmAtomName;
    double charmmPartCharge;
  };

  struct param19 {
    double charmmMass;
    double charmmMinEnergy;
    double charmmMinRadius;
  };

  // set atom.charmmMass, atom.charmmMinEnergy, atom.charmmMinRadius
  static void setCharmmParams(cAtom *atom, cAminoResidue::eAtomLabel atomLabel);

 private:
  static const std::map<cAminoResidue::eAtomLabel, toph19> mCharmmAtomsALA;
  static const std::map<cAminoResidue::eAtomLabel, toph19> mCharmmAtomsARG;
  static const std::map<cAminoResidue::eAtomLabel, toph19> mCharmmAtomsASN;
  static const std::map<cAminoResidue::eAtomLabel, toph19> mCharmmAtomsASP;
  static const std::map<cAminoResidue::eAtomLabel, toph19> mCharmmAtomsCYS;
  static const std::map<cAminoResidue::eAtomLabel, toph19> mCharmmAtomsGLN;
  static const std::map<cAminoResidue::eAtomLabel, toph19> mCharmmAtomsGLU;
  static const std::map<cAminoResidue::eAtomLabel, toph19> mCharmmAtomsGLY;
  static const std::map<cAminoResidue::eAtomLabel, toph19> mCharmmAtomsHIS;
  static const std::map<cAminoResidue::eAtomLabel, toph19> mCharmmAtomsILE;
  static const std::map<cAminoResidue::eAtomLabel, toph19> mCharmmAtomsLEU;
  static const std::map<cAminoResidue::eAtomLabel, toph19> mCharmmAtomsLYS;
  static const std::map<cAminoResidue::eAtomLabel, toph19> mCharmmAtomsMET;
  static const std::map<cAminoResidue::eAtomLabel, toph19> mCharmmAtomsPHE;
  static const std::map<cAminoResidue::eAtomLabel, toph19> mCharmmAtomsPRO;
  static const std::map<cAminoResidue::eAtomLabel, toph19> mCharmmAtomsSER;
  static const std::map<cAminoResidue::eAtomLabel, toph19> mCharmmAtomsTHR;
  static const std::map<cAminoResidue::eAtomLabel, toph19> mCharmmAtomsTRP;
  static const std::map<cAminoResidue::eAtomLabel, toph19> mCharmmAtomsTYR;
  static const std::map<cAminoResidue::eAtomLabel, toph19> mCharmmAtomsVAL;
  static const std::map<cAminoResidue::eAtomLabel, toph19> mCharmmAtomsSEC;
  static const std::map<cAminoResidue::eAtomLabel, toph19> mCharmmAtomsMSE;

  static const std::map<eAtomName, param19> mAtomParams;
};
