#include "cCharmmParams.hpp"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>


const std::map<cAminoResidue::eAtomLabel, cCharmmParams::toph19> cCharmmParams::mCharmmAtomsALA = {
  {cAminoResidue::eAtomLabel::N, {NH1, -0.35}},
  {cAminoResidue::eAtomLabel::H, {H, 0.25}},
  {cAminoResidue::eAtomLabel::CA, {CH1E, 0.10}},

  {cAminoResidue::eAtomLabel::CB, {CH3E, 0.00}},

  {cAminoResidue::eAtomLabel::C, {C, 0.55}},
  {cAminoResidue::eAtomLabel::O, {O, -0.55}},

  {cAminoResidue::eAtomLabel::HA, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HB3, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HB2, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HB1, {H, 0.00}},
  {cAminoResidue::eAtomLabel::OXT, {O, 0.00}}, // ??? check this
  {cAminoResidue::eAtomLabel::HT1, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HT2, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HT3, {H, 0.00}},
};

const std::map<cAminoResidue::eAtomLabel, cCharmmParams::toph19> cCharmmParams::mCharmmAtomsARG = {
  {cAminoResidue::eAtomLabel::N, {NH1, -0.35}},
  {cAminoResidue::eAtomLabel::H, {H, 0.25}},
  {cAminoResidue::eAtomLabel::CA, {CH1E, 0.10}},

  {cAminoResidue::eAtomLabel::CB, {CH2E, 0.00}},
  {cAminoResidue::eAtomLabel::CG, {CH2E, 0.00}},

  {cAminoResidue::eAtomLabel::CD, {CH2E, 0.10}},
  {cAminoResidue::eAtomLabel::NE, {NH1, -0.40}},
  {cAminoResidue::eAtomLabel::HE, {H, 0.30}},
  {cAminoResidue::eAtomLabel::CZ, {C, 0.50}},

  {cAminoResidue::eAtomLabel::NH1, {NC2, -0.45}},
  {cAminoResidue::eAtomLabel::HH11, {HC, 0.35}},
  {cAminoResidue::eAtomLabel::HH12, {HC, 0.35}},

  {cAminoResidue::eAtomLabel::NH2, {NC2, -0.45}},
  {cAminoResidue::eAtomLabel::HH21, {HC, 0.35}},
  {cAminoResidue::eAtomLabel::HH22, {HC, 0.35}},

  {cAminoResidue::eAtomLabel::C, {C, 0.55}},
  {cAminoResidue::eAtomLabel::O, {O, -0.55}},

  {cAminoResidue::eAtomLabel::HA, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HB2, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HB1, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HD2, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HD1, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HG2, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HG1, {H, 0.00}},
  {cAminoResidue::eAtomLabel::OXT, {O, 0.00}}, // ??? check this
  {cAminoResidue::eAtomLabel::HT1, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HT2, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HT3, {H, 0.00}},
};

const std::map<cAminoResidue::eAtomLabel, cCharmmParams::toph19> cCharmmParams::mCharmmAtomsASN = {
  {cAminoResidue::eAtomLabel::N, {NH1, -0.35}},
  {cAminoResidue::eAtomLabel::H, {H, 0.25}},
  {cAminoResidue::eAtomLabel::CA, {CH1E, 0.10}},

  {cAminoResidue::eAtomLabel::CB, {CH2E, 0.00}},

  {cAminoResidue::eAtomLabel::CG, {C, 0.55}},
  {cAminoResidue::eAtomLabel::OD1, {O, -0.55}},

  {cAminoResidue::eAtomLabel::ND2, {NH2, -0.60}},
  {cAminoResidue::eAtomLabel::HD21, {H, 0.30}},
  {cAminoResidue::eAtomLabel::HD22, {H, 0.30}},

  {cAminoResidue::eAtomLabel::C, {C, 0.55}},
  {cAminoResidue::eAtomLabel::O, {O, -0.55}},

  {cAminoResidue::eAtomLabel::HA, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HB2, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HB1, {H, 0.00}},
  {cAminoResidue::eAtomLabel::OXT, {O, 0.00}}, // ??? check this
  {cAminoResidue::eAtomLabel::HT1, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HT2, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HT3, {H, 0.00}},
};

const std::map<cAminoResidue::eAtomLabel, cCharmmParams::toph19> cCharmmParams::mCharmmAtomsASP = {
  {cAminoResidue::eAtomLabel::N, {NH1, -0.35}},
  {cAminoResidue::eAtomLabel::H, {H, 0.25}},
  {cAminoResidue::eAtomLabel::CA, {CH1E, 0.10}},

  {cAminoResidue::eAtomLabel::CB, {CH2E, -0.16}},
  {cAminoResidue::eAtomLabel::CG, {C, 0.36}},
  {cAminoResidue::eAtomLabel::OD1, {OC, -0.60}},
  {cAminoResidue::eAtomLabel::OD2, {OC, -0.60}},

  {cAminoResidue::eAtomLabel::C, {C, 0.55}},
  {cAminoResidue::eAtomLabel::O, {O, -0.55}},

  {cAminoResidue::eAtomLabel::HA, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HB2, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HB1, {H, 0.00}},
  {cAminoResidue::eAtomLabel::OXT, {O, 0.00}}, // ??? check this
  {cAminoResidue::eAtomLabel::HT1, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HT2, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HT3, {H, 0.00}},
};

const std::map<cAminoResidue::eAtomLabel, cCharmmParams::toph19> cCharmmParams::mCharmmAtomsCYS = {
  {cAminoResidue::eAtomLabel::N, {NH1, -0.35}},
  {cAminoResidue::eAtomLabel::H, {H, 0.25}},
  {cAminoResidue::eAtomLabel::CA, {CH1E, 0.10}},

  {cAminoResidue::eAtomLabel::CB, {CH2E, 0.19}},
  {cAminoResidue::eAtomLabel::SG, {SH1E, -0.19}},

  {cAminoResidue::eAtomLabel::C, {C, 0.55}},
  {cAminoResidue::eAtomLabel::O, {O, -0.55}},

  {cAminoResidue::eAtomLabel::HA, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HB2, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HB1, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HG, {H, 0.00}},
  {cAminoResidue::eAtomLabel::OXT, {O, 0.00}}, // ??? check this
  {cAminoResidue::eAtomLabel::HT1, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HT2, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HT3, {H, 0.00}},
};

// FIXME: This is not from charmm but just a copy of mCharmmAtomsCYS
const std::map<cAminoResidue::eAtomLabel, cCharmmParams::toph19> cCharmmParams::mCharmmAtomsSEC = {
  {cAminoResidue::eAtomLabel::N, {NH1, -0.35}},
  {cAminoResidue::eAtomLabel::H, {H, 0.25}},
  {cAminoResidue::eAtomLabel::CA, {CH1E, 0.10}},

  {cAminoResidue::eAtomLabel::CB, {CH2E, 0.19}},
  {cAminoResidue::eAtomLabel::SE, {SE, -0.19}},

  {cAminoResidue::eAtomLabel::C, {C, 0.55}},
  {cAminoResidue::eAtomLabel::O, {O, -0.55}},

  {cAminoResidue::eAtomLabel::HA, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HB2, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HB1, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HG, {H, 0.00}},
  {cAminoResidue::eAtomLabel::OXT, {O, 0.00}}, // ??? check this
  {cAminoResidue::eAtomLabel::HT1, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HT2, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HT3, {H, 0.00}},
};

const std::map<cAminoResidue::eAtomLabel, cCharmmParams::toph19> cCharmmParams::mCharmmAtomsGLN = {
  {cAminoResidue::eAtomLabel::N, {NH1, -0.35}},
  {cAminoResidue::eAtomLabel::H, {H, 0.25}},
  {cAminoResidue::eAtomLabel::CA, {CH1E, 0.10}},

  {cAminoResidue::eAtomLabel::CB, {CH2E, 0.00}},
  {cAminoResidue::eAtomLabel::CG, {CH2E, 0.00}},

  {cAminoResidue::eAtomLabel::CD, {C, 0.55}},
  {cAminoResidue::eAtomLabel::OE1, {O, -0.55}},

  {cAminoResidue::eAtomLabel::NE2, {NH2, -0.60}},
  {cAminoResidue::eAtomLabel::HE21, {H, 0.30}},
  {cAminoResidue::eAtomLabel::HE22, {H, 0.30}},

  {cAminoResidue::eAtomLabel::C, {C, 0.55}},
  {cAminoResidue::eAtomLabel::O, {O, -0.55}},

  {cAminoResidue::eAtomLabel::HA, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HB2, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HB1, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HG2, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HG1, {H, 0.00}},
  {cAminoResidue::eAtomLabel::OXT, {O, 0.00}}, // ??? check this
  {cAminoResidue::eAtomLabel::HT1, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HT2, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HT3, {H, 0.00}},
};

const std::map<cAminoResidue::eAtomLabel, cCharmmParams::toph19> cCharmmParams::mCharmmAtomsGLU = {
  {cAminoResidue::eAtomLabel::N, {NH1, -0.35}},
  {cAminoResidue::eAtomLabel::H, {H, 0.25}},
  {cAminoResidue::eAtomLabel::CA, {CH1E, 0.10}},

  {cAminoResidue::eAtomLabel::CB, {CH2E, 0.00}},

  {cAminoResidue::eAtomLabel::CG, {CH2E, -0.16}},
  {cAminoResidue::eAtomLabel::CD, {C, 0.36}},
  {cAminoResidue::eAtomLabel::OE1, {OC, -0.60}},
  {cAminoResidue::eAtomLabel::OE2, {OC, -0.60}},

  {cAminoResidue::eAtomLabel::C, {C, 0.55}},
  {cAminoResidue::eAtomLabel::O, {O, -0.55}},

  {cAminoResidue::eAtomLabel::HA, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HB2, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HB1, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HG2, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HG1, {H, 0.00}},
  {cAminoResidue::eAtomLabel::OXT, {O, 0.00}}, // ??? check this
  {cAminoResidue::eAtomLabel::HT1, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HT2, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HT3, {H, 0.00}},
};

const std::map<cAminoResidue::eAtomLabel, cCharmmParams::toph19> cCharmmParams::mCharmmAtomsGLY = {
  {cAminoResidue::eAtomLabel::N, {NH1, -0.35}},
  {cAminoResidue::eAtomLabel::H, {H, 0.25}},
  {cAminoResidue::eAtomLabel::CA, {CH2E, 0.10}},

  {cAminoResidue::eAtomLabel::C, {C, 0.55}},
  {cAminoResidue::eAtomLabel::O, {O, -0.55}},

  {cAminoResidue::eAtomLabel::HA2, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HA1, {H, 0.00}},
  {cAminoResidue::eAtomLabel::OXT, {O, 0.00}}, // ??? check this
  {cAminoResidue::eAtomLabel::HT1, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HT2, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HT3, {H, 0.00}},
};

const std::map<cAminoResidue::eAtomLabel, cCharmmParams::toph19> cCharmmParams::mCharmmAtomsHIS = {
  {cAminoResidue::eAtomLabel::N, {NH1, -0.35}},
  {cAminoResidue::eAtomLabel::H, {H, 0.25}},
  {cAminoResidue::eAtomLabel::CA, {CH1E, 0.10}},

  {cAminoResidue::eAtomLabel::CB, {CH2E, 0.00}},

  {cAminoResidue::eAtomLabel::CG, {C, 0.10}},
  {cAminoResidue::eAtomLabel::ND1, {NH1, -0.40}},
  {cAminoResidue::eAtomLabel::HD1, {H, 0.30}},

  {cAminoResidue::eAtomLabel::CD2, {CR1E, 0.10}},
  {cAminoResidue::eAtomLabel::NE2, {NR, -0.40}},
  {cAminoResidue::eAtomLabel::CE1, {CR1E, 0.30}},

  {cAminoResidue::eAtomLabel::C, {C, 0.55}},
  {cAminoResidue::eAtomLabel::O, {O, -0.55}},

  {cAminoResidue::eAtomLabel::HA, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HB2, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HB1, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HD2, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HE1, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HE2, {H, 0.00}},
  {cAminoResidue::eAtomLabel::OXT, {O, 0.00}}, // ??? check this
  {cAminoResidue::eAtomLabel::HT1, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HT2, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HT3, {H, 0.00}},
};

const std::map<cAminoResidue::eAtomLabel, cCharmmParams::toph19> cCharmmParams::mCharmmAtomsILE = {
  {cAminoResidue::eAtomLabel::N, {NH1, -0.35}},
  {cAminoResidue::eAtomLabel::H, {H, 0.25}},
  {cAminoResidue::eAtomLabel::CA, {CH1E, 0.10}},

  {cAminoResidue::eAtomLabel::CB, {CH1E, 0.00}},
  {cAminoResidue::eAtomLabel::CG2, {CH3E, 0.00}},

  {cAminoResidue::eAtomLabel::CG1, {CH2E, 0.00}},
  {cAminoResidue::eAtomLabel::CD1, {CH3E, 0.00}},

  {cAminoResidue::eAtomLabel::C, {C, 0.55}},
  {cAminoResidue::eAtomLabel::O, {O, -0.55}},

  {cAminoResidue::eAtomLabel::HA, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HB, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HD13, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HD12, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HD11, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HG12, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HG23, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HG11, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HG22, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HG21, {H, 0.00}},
  {cAminoResidue::eAtomLabel::OXT, {O, 0.00}}, // ??? check this
  {cAminoResidue::eAtomLabel::HT1, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HT2, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HT3, {H, 0.00}},
};

const std::map<cAminoResidue::eAtomLabel, cCharmmParams::toph19> cCharmmParams::mCharmmAtomsLEU = {
  {cAminoResidue::eAtomLabel::N, {NH1, -0.35}},
  {cAminoResidue::eAtomLabel::H, {H, 0.25}},
  {cAminoResidue::eAtomLabel::CA, {CH1E, 0.10}},

  {cAminoResidue::eAtomLabel::CB, {CH2E, 0.00}},
  {cAminoResidue::eAtomLabel::CG, {CH1E, 0.00}},
  {cAminoResidue::eAtomLabel::CD1, {CH3E, 0.00}},
  {cAminoResidue::eAtomLabel::CD2, {CH3E, 0.00}},

  {cAminoResidue::eAtomLabel::C, {C, 0.55}},
  {cAminoResidue::eAtomLabel::O, {O, -0.55}},

  {cAminoResidue::eAtomLabel::HA, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HB2, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HB1, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HD13, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HD23, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HD12, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HD22, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HD11, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HD21, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HG, {H, 0.00}},
  {cAminoResidue::eAtomLabel::OXT, {O, 0.00}}, // ??? check this,
  {cAminoResidue::eAtomLabel::HT1, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HT2, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HT3, {H, 0.00}},
};

const std::map<cAminoResidue::eAtomLabel, cCharmmParams::toph19> cCharmmParams::mCharmmAtomsLYS = {
  {cAminoResidue::eAtomLabel::N, {NH1, -0.35}},
  {cAminoResidue::eAtomLabel::H, {H, 0.25}},
  {cAminoResidue::eAtomLabel::CA, {CH1E, 0.10}},

  {cAminoResidue::eAtomLabel::CB, {CH2E, 0.00}},
  {cAminoResidue::eAtomLabel::CG, {CH2E, 0.00}},
  {cAminoResidue::eAtomLabel::CD, {CH2E, 0.00}},

  {cAminoResidue::eAtomLabel::CE, {CH2E, 0.25}},
  {cAminoResidue::eAtomLabel::NZ, {NH3, -0.30}},
  {cAminoResidue::eAtomLabel::HZ1, {HC, 0.35}},
  {cAminoResidue::eAtomLabel::HZ2, {HC, 0.35}},
  {cAminoResidue::eAtomLabel::HZ3, {HC, 0.35}},

  {cAminoResidue::eAtomLabel::C, {C, 0.55}},
  {cAminoResidue::eAtomLabel::O, {O, -0.55}},

  {cAminoResidue::eAtomLabel::HA, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HB2, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HB1, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HD2, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HD1, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HE2, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HE1, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HG2, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HG1, {H, 0.00}},
  {cAminoResidue::eAtomLabel::OXT, {O, 0.00}}, // ??? check this
  {cAminoResidue::eAtomLabel::HT1, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HT2, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HT3, {H, 0.00}},
};

const std::map<cAminoResidue::eAtomLabel, cCharmmParams::toph19> cCharmmParams::mCharmmAtomsMET = {
  {cAminoResidue::eAtomLabel::N, {NH1, -0.35}},
  {cAminoResidue::eAtomLabel::H, {H, 0.25}},
  {cAminoResidue::eAtomLabel::CA, {CH1E, 0.10}},

  {cAminoResidue::eAtomLabel::CB, {CH2E, 0.00}},

  {cAminoResidue::eAtomLabel::CG, {CH2E, 0.06}},
  {cAminoResidue::eAtomLabel::SD, {S, -0.12}},
  {cAminoResidue::eAtomLabel::CE, {CH3E, 0.06}},

  {cAminoResidue::eAtomLabel::C, {C, 0.55}},
  {cAminoResidue::eAtomLabel::O, {O, -0.55}},

  {cAminoResidue::eAtomLabel::HA, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HB2, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HB1, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HE3, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HE2, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HE1, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HG2, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HG1, {H, 0.00}},
  {cAminoResidue::eAtomLabel::OXT, {O, 0.00}}, // ??? check this
  {cAminoResidue::eAtomLabel::HT1, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HT2, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HT3, {H, 0.00}},
};

// FIXME: This is not from charmm but just a copy of mCharmmAtomsMET
const std::map<cAminoResidue::eAtomLabel, cCharmmParams::toph19> cCharmmParams::mCharmmAtomsMSE = {
  {cAminoResidue::eAtomLabel::N, {NH1, -0.35}},
  {cAminoResidue::eAtomLabel::H, {H, 0.25}},
  {cAminoResidue::eAtomLabel::CA, {CH1E, 0.10}},

  {cAminoResidue::eAtomLabel::CB, {CH2E, 0.00}},

  {cAminoResidue::eAtomLabel::CG, {CH2E, 0.06}},
  {cAminoResidue::eAtomLabel::SE, {SE, -0.12}},
  {cAminoResidue::eAtomLabel::CE, {CH3E, 0.06}},

  {cAminoResidue::eAtomLabel::C, {C, 0.55}},
  {cAminoResidue::eAtomLabel::O, {O, -0.55}},

  {cAminoResidue::eAtomLabel::HA, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HB2, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HB1, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HE3, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HE2, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HE1, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HG2, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HG1, {H, 0.00}},
  {cAminoResidue::eAtomLabel::OXT, {O, 0.00}}, // ??? check this
  {cAminoResidue::eAtomLabel::HT1, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HT2, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HT3, {H, 0.00}},
};

const std::map<cAminoResidue::eAtomLabel, cCharmmParams::toph19> cCharmmParams::mCharmmAtomsPHE = {
  {cAminoResidue::eAtomLabel::N, {NH1, -0.35}},
  {cAminoResidue::eAtomLabel::H, {H, 0.25}},
  {cAminoResidue::eAtomLabel::CA, {CH1E, 0.10}},

  {cAminoResidue::eAtomLabel::CB, {CH2E, 0.00}},
  {cAminoResidue::eAtomLabel::CG, {C, 0.00}},
  {cAminoResidue::eAtomLabel::CD1, {CR1E, 0.00}},
  {cAminoResidue::eAtomLabel::CD2, {CR1E, 0.00}},

  {cAminoResidue::eAtomLabel::CE1, {CR1E, 0.00}},
  {cAminoResidue::eAtomLabel::CE2, {CR1E, 0.00}},
  {cAminoResidue::eAtomLabel::CZ, {CR1E, 0.00}},

  {cAminoResidue::eAtomLabel::C, {C, 0.55}},
  {cAminoResidue::eAtomLabel::O, {O, -0.55}},

  {cAminoResidue::eAtomLabel::HA, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HB2, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HB1, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HD1, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HD2, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HE1, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HE2, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HZ, {H, 0.00}},
  {cAminoResidue::eAtomLabel::OXT, {O, 0.00}}, // ??? check this
  {cAminoResidue::eAtomLabel::HT1, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HT2, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HT3, {H, 0.00}},
};

const std::map<cAminoResidue::eAtomLabel, cCharmmParams::toph19> cCharmmParams::mCharmmAtomsPRO = {
  {cAminoResidue::eAtomLabel::N, {N, -0.20}},
  {cAminoResidue::eAtomLabel::CD, {CH2E, 0.10}},
  {cAminoResidue::eAtomLabel::CA, {CH1E, 0.10}},

  {cAminoResidue::eAtomLabel::CB, {CH2E, 0.00}},
  {cAminoResidue::eAtomLabel::CG, {CH2E, 0.00}},

  {cAminoResidue::eAtomLabel::C, {C, 0.55}},
  {cAminoResidue::eAtomLabel::O, {O, -0.55}},

  {cAminoResidue::eAtomLabel::HA, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HB2, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HB1, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HG2, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HG1, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HD2, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HD1, {H, 0.00}},
  {cAminoResidue::eAtomLabel::OXT, {O, 0.00}}, // ??? check this
  {cAminoResidue::eAtomLabel::HT1, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HT2, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HT3, {H, 0.00}},
};

const std::map<cAminoResidue::eAtomLabel, cCharmmParams::toph19> cCharmmParams::mCharmmAtomsSER = {
  {cAminoResidue::eAtomLabel::N, {NH1, -0.35}},
  {cAminoResidue::eAtomLabel::H, {H, 0.25}},
  {cAminoResidue::eAtomLabel::CA, {CH1E, 0.10}},

  {cAminoResidue::eAtomLabel::CB, {CH2E, 0.25}},
  {cAminoResidue::eAtomLabel::OG, {OH1, -0.65}},
  {cAminoResidue::eAtomLabel::HG, {H, 0.40}},

  {cAminoResidue::eAtomLabel::C, {C, 0.55}},
  {cAminoResidue::eAtomLabel::O, {O, -0.55}},

  {cAminoResidue::eAtomLabel::HA, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HB2, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HB1, {H, 0.00}},
  {cAminoResidue::eAtomLabel::OXT, {O, 0.00}}, // ??? check this
  {cAminoResidue::eAtomLabel::HT1, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HT2, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HT3, {H, 0.00}},
};

const std::map<cAminoResidue::eAtomLabel, cCharmmParams::toph19> cCharmmParams::mCharmmAtomsTHR = {
  {cAminoResidue::eAtomLabel::N, {NH1, -0.35}},
  {cAminoResidue::eAtomLabel::H, {H, 0.25}},
  {cAminoResidue::eAtomLabel::CA, {CH1E, 0.10}},

  {cAminoResidue::eAtomLabel::CB, {CH1E, 0.25}},
  {cAminoResidue::eAtomLabel::OG1, {OH1, -0.65}},
  {cAminoResidue::eAtomLabel::HG1, {H, 0.40}},

  {cAminoResidue::eAtomLabel::CG2, {CH3E, 0.00}},

  {cAminoResidue::eAtomLabel::C, {C, 0.55}},
  {cAminoResidue::eAtomLabel::O, {O, -0.55}},

  {cAminoResidue::eAtomLabel::HA, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HB, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HG23, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HG22, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HG21, {H, 0.00}},
  {cAminoResidue::eAtomLabel::OXT, {O, 0.00}}, // ??? check this
  {cAminoResidue::eAtomLabel::HT1, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HT2, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HT3, {H, 0.00}},
};

const std::map<cAminoResidue::eAtomLabel, cCharmmParams::toph19> cCharmmParams::mCharmmAtomsTRP = {
  {cAminoResidue::eAtomLabel::N, {NH1, -0.35}},
  {cAminoResidue::eAtomLabel::H, {H, 0.25}},
  {cAminoResidue::eAtomLabel::CA, {CH1E, 0.10}},

  {cAminoResidue::eAtomLabel::CB, {CH2E, 0.00}},

  {cAminoResidue::eAtomLabel::CG, {C, -0.03}},
  {cAminoResidue::eAtomLabel::CD2, {C, 0.10}},
  {cAminoResidue::eAtomLabel::CE2, {C, -0.04}},
  {cAminoResidue::eAtomLabel::CE3, {CR1E, -0.03}},

  {cAminoResidue::eAtomLabel::CD1, {CR1E, 0.06}},
  {cAminoResidue::eAtomLabel::NE1, {NH1, -0.36}},
  {cAminoResidue::eAtomLabel::HE1, {H, 0.30}},

  {cAminoResidue::eAtomLabel::CZ2, {CR1E, 0.00}},
  {cAminoResidue::eAtomLabel::CZ3, {CR1E, 0.00}},
  {cAminoResidue::eAtomLabel::CH2, {CR1E, 0.00}},

  {cAminoResidue::eAtomLabel::C, {C, 0.55}},
  {cAminoResidue::eAtomLabel::O, {O, -0.55}},

  {cAminoResidue::eAtomLabel::HA, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HB2, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HB1, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HD1, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HE3, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HH2, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HZ2, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HZ3, {H, 0.00}},
  {cAminoResidue::eAtomLabel::OXT, {O, 0.00}}, // ??? check this
  {cAminoResidue::eAtomLabel::HT1, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HT2, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HT3, {H, 0.00}},
};

const std::map<cAminoResidue::eAtomLabel, cCharmmParams::toph19> cCharmmParams::mCharmmAtomsTYR = {
  {cAminoResidue::eAtomLabel::N, {NH1, -0.35}},
  {cAminoResidue::eAtomLabel::H, {H, 0.25}},
  {cAminoResidue::eAtomLabel::CA, {CH1E, 0.10}},

  {cAminoResidue::eAtomLabel::CB, {CH2E, 0.00}},
  {cAminoResidue::eAtomLabel::CG, {C, 0.00}},

  {cAminoResidue::eAtomLabel::CD1, {CR1E, 0.00}},
  {cAminoResidue::eAtomLabel::CE1, {CR1E, 0.00}},

  {cAminoResidue::eAtomLabel::CD2, {CR1E, 0.00}},
  {cAminoResidue::eAtomLabel::CE2, {CR1E, 0.00}},

  {cAminoResidue::eAtomLabel::CZ, {C, 0.25}},
  {cAminoResidue::eAtomLabel::OH, {OH1, -0.65}},
  {cAminoResidue::eAtomLabel::HH, {H, 0.40}},

  {cAminoResidue::eAtomLabel::C, {C, 0.55}},
  {cAminoResidue::eAtomLabel::O, {O, -0.55}},

  {cAminoResidue::eAtomLabel::HA, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HB2, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HB1, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HD1, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HD2, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HE1, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HE2, {H, 0.00}},
  {cAminoResidue::eAtomLabel::OXT, {O, 0.00}}, // ??? check this
  {cAminoResidue::eAtomLabel::HT1, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HT2, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HT3, {H, 0.00}},
};

const std::map<cAminoResidue::eAtomLabel, cCharmmParams::toph19> cCharmmParams::mCharmmAtomsVAL = {
  {cAminoResidue::eAtomLabel::N, {NH1, -0.35}},
  {cAminoResidue::eAtomLabel::H, {H, 0.25}},
  {cAminoResidue::eAtomLabel::CA, {CH1E, 0.10}},

  {cAminoResidue::eAtomLabel::CB, {CH1E, 0.00}},
  {cAminoResidue::eAtomLabel::CG1, {CH3E, 0.00}},
  {cAminoResidue::eAtomLabel::CG2, {CH3E, 0.00}},

  {cAminoResidue::eAtomLabel::C, {C, 0.55}},
  {cAminoResidue::eAtomLabel::O, {O, -0.55}},

  {cAminoResidue::eAtomLabel::HA, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HB, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HG13, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HG23, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HG12, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HG22, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HG11, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HG21, {H, 0.00}},
  {cAminoResidue::eAtomLabel::OXT, {O, 0.00}}, // ??? check this
  {cAminoResidue::eAtomLabel::HT1, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HT2, {H, 0.00}},
  {cAminoResidue::eAtomLabel::HT3, {H, 0.00}},
};

const std::map<cCharmmParams::eAtomName, cCharmmParams::param19> cCharmmParams::mAtomParams = {
  {H, {1.00800, -0.0498, 0.8000}},
  {HC, {1.00800, -0.0498, 0.6000}},
  {HA, {1.00800, -0.0450 , 1.4680}},
  {HT, {1.00800, -0.0498, 0.8000}},

  {C, {12.01100, -0.1200, 2.100}},
  {CH1E, {13.01900, -0.0486, 2.365}},
  {CH2E, {14.02700, -0.1142, 2.235}},
  {CH3E, {15.03500, -0.1811, 2.165}},
  {CR1E, {13.01900, -0.1200, 2.100}},
  {CT, {12.01100, -0.0262, 2.490}},
  {CM, {12.01100, -0.0262, 2.490}},

  {N, {14.00670, -0.2384, 1.6000}},
  {NR, {14.00670, -0.2384, 1.6000}},
  {NP, {14.00670, -0.2384, 1.6000}},
  {NH1E, {15.01470, -0.2384, 1.6000}}, // ??? check this
  {NH2E, {16.02270, -0.2384, 1.6000}}, // ??? check this
  {NH3E, {17.03070, -0.2384, 1.6000}}, // ??? check this
  {NC2E, {16.02270, -0.2384, 1.6000}}, // ??? check this
  {NH1, {14.00670, -0.2384, 1.6000}},
  {NH2, {14.00670, -0.2384, 1.6000}},
  {NH3, {14.00670, -0.2384, 1.6000}},
  {NC2, {14.00670, -0.2384, 1.6000}},

  {O, {15.99940, -0.1591, 1.6000}},
  {OC, {15.99940, -0.6469, 1.6000}},
  {OH1E, {17.00740, -0.1591, 1.6000}}, // ??? check this
  {OH2E, {18.01540, -0.1591, 1.6000}}, // ??? check this
  {OH1, {15.99940, -0.1591, 1.6000}},
  {OH2, {15.99940, -0.0758, 1.7398}},
  {OM, {15.99940, -0.1591, 1.6000}},
  {OT, {15.99940, -0.1591, 1.6000}},
  {OS, {15.99940, -0.1591, 1.6000}},

  {LP, {0.0, -0.04598, 0.2245}},
  {FE, {55.84700, 0.000, 0.6500}},
  {S, {32.06000, -0.0430, 1.890}},
  {SH1E, {33.06800, -0.0430, 1.890}},
  {SE, {33.06800, -0.0430, 1.890}}, // FIXME: This is not from charmm

  {CAL, {40.08000, -0.120000, 1.710000}}
};

// set atom.CharmmMass, atom.CharmmMinEnergy, atom.CharmmMinRadius
void cCharmmParams::setCharmmParams(cAtom *atom, cAminoResidue::eAtomLabel atomLabel) {
  cAminoResidue *residue
    = dynamic_cast<cAminoResidue *>(atom->residue);

  const std::map<cAminoResidue::eAtomLabel, toph19> *mCharmmAtoms;
  switch (residue->aminoType()) {
    case cAminoResidue::ALA:
      mCharmmAtoms = &mCharmmAtomsALA;
      break;
    case cAminoResidue::ARG:
      mCharmmAtoms = &mCharmmAtomsARG;
      break;
    case cAminoResidue::ASN:
      mCharmmAtoms = &mCharmmAtomsASN;
      break;
    case cAminoResidue::ASP:
      mCharmmAtoms = &mCharmmAtomsASP;
      break;
    case cAminoResidue::CYS:
      mCharmmAtoms = &mCharmmAtomsCYS;
      break;
    case cAminoResidue::SEC:
      mCharmmAtoms = &mCharmmAtomsSEC;
      break;
    case cAminoResidue::GLN:
      mCharmmAtoms = &mCharmmAtomsGLN;
      break;
    case cAminoResidue::GLU:
      mCharmmAtoms = &mCharmmAtomsGLU;
      break;
    case cAminoResidue::GLY:
      mCharmmAtoms = &mCharmmAtomsGLY;
      break;
    case cAminoResidue::HIS:
      mCharmmAtoms = &mCharmmAtomsHIS;
      break;
    case cAminoResidue::ILE:
      mCharmmAtoms = &mCharmmAtomsILE;
      break;
    case cAminoResidue::LEU:
      mCharmmAtoms = &mCharmmAtomsLEU;
      break;
    case cAminoResidue::LYS:
      mCharmmAtoms = &mCharmmAtomsLYS;
      break;
    case cAminoResidue::MET:
      mCharmmAtoms = &mCharmmAtomsMET;
      break;
    case cAminoResidue::MSE:
      mCharmmAtoms = &mCharmmAtomsMSE;
      break;
    case cAminoResidue::PHE:
      mCharmmAtoms = &mCharmmAtomsPHE;
      break;
    case cAminoResidue::PRO:
      mCharmmAtoms = &mCharmmAtomsPRO;
      break;
    case cAminoResidue::SER:
      mCharmmAtoms = &mCharmmAtomsSER;
      break;
    case cAminoResidue::THR:
      mCharmmAtoms = &mCharmmAtomsTHR;
      break;
    case cAminoResidue::TRP:
      mCharmmAtoms = &mCharmmAtomsTRP;
      break;
    case cAminoResidue::TYR:
      mCharmmAtoms = &mCharmmAtomsTYR;
      break;
    case cAminoResidue::VAL:
      mCharmmAtoms = &mCharmmAtomsVAL;
      break;
    default:
      fprintf(stderr, "Unknown type of amino residue: %s\n",
                                        residue->resName.c_str());
      exit(1);
  }
  std::map<cAminoResidue::eAtomLabel, toph19>::const_iterator it
    = mCharmmAtoms->find(atomLabel);

  assert(it != mCharmmAtoms->end());
  const toph19 &toph19_ = it->second;
  if (toph19_.charmmAtomName == DUMMY)
    return;

  atom->setCharmmPartCharge(toph19_.charmmPartCharge);
  std::map<eAtomName, param19>::const_iterator it_
    = mAtomParams.find(toph19_.charmmAtomName);

  assert(it_ != mAtomParams.end());
  const param19 &param19_ = it_->second;
  atom->setCharmmMass(param19_.charmmMass);
  atom->setCharmmMinEnergy(param19_.charmmMinEnergy);
  atom->setCharmmMinRadius(param19_.charmmMinRadius);
}
