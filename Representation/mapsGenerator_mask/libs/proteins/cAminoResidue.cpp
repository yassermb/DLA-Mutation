#include "cAminoResidue.hpp"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <algorithm>
#include <stdexcept>

#include "cRotamerLibrary.hpp"
#include "cCharmmParams.hpp"
#include "mathFunctions.hpp"


const cAminoResidue::residueRebuildingContext cAminoResidue::ALAresidueRebuildingData[] = {
  {{N, C, CA, HA}, 1.090000, 109.000000, -121.435662, 0},
  {{N, C, CA, CB}, 1.521000, 110.500000, +122.968111, 0},

  {{N, CA, CB, HB1}, 1.090000, 109.000000, 180.0, 0},
  {{N, CA, CB, HB2}, 1.090000, 109.000000, 300.0, 0},
  {{N, CA, CB, HB3}, 1.090000, 109.000000, 60.0, 0},
};

const cAminoResidue::residueRebuildingContext cAminoResidue::ARGresidueRebuildingData[] = {
  {{N, C, CA, HA}, 1.090000, 109.000000, -121.435662, 0},
  {{N, C, CA, CB}, 1.530000, 110.100000, +122.815261, 0},

  {{N, CA, CB, CG}, 1.520000, 114.100000, 0, 1},
  {{CG, CA, CB, HB1}, 1.090000, 109.000000, +120.800876, 0},
  {{CG, CA, CB, HB2}, 1.090000, 109.000000, -120.800876, 0},

  {{CA, CB, CG, CD}, 1.520000, 111.300000, 0, 2},
  {{CD, CB, CG, HG1}, 1.090000, 109.000000, +119.014578, 0},
  {{CD, CB, CG, HG2}, 1.090000, 109.000000, -119.014578, 0},

  {{CB, CG, CD, NE}, 1.460000, 112.000000, 0, 3},
  {{NE, CG, CD, HD1}, 1.090000, 109.000000, +119.446266, 0},
  {{NE, CG, CD, HD2}, 1.090000, 109.000000, -119.446266, 0},

  {{CG, CD, NE, CZ}, 1.329000, 124.200000, 0, 4},
  {{CZ, CD, NE, HE}, 1.020000, 117.900000, 179.999999, 0},

  {{CD, NE, CZ, NH1}, 1.326000, 120.000000, 180.0, 0},
  {{NE, CZ, NH1, HH11}, 1.020000, 120.000000, 0.0, 0},
  {{NE, CZ, NH1, HH12}, 1.020000, 120.000000, 180.0, 0},

  {{CD, NE, CZ, NH2}, 1.326000, 120.000000, 0.0, 0},
  {{NE, CZ, NH2, HH21}, 1.020000, 120.000000, 0.0, 0},
  {{NE, CZ, NH2, HH22}, 1.020000, 120.000000, 180.0, 0},
};

const cAminoResidue::residueRebuildingContext cAminoResidue::ASNresidueRebuildingData[] = {
  {{N, C, CA, HA}, 1.090000, 109.000000, -122.087018, 0},
  {{N, C, CA, CB}, 1.530000, 110.100000, +123.512487, 0},

  {{N, CA, CB, CG}, 1.516000, 112.600000, 0, 1},
  {{CG, CA, CB, HB1}, 1.090000, 109.000000, +119.823949, 0},
  {{CG, CA, CB, HB2}, 1.090000, 109.000000, -119.823949, 0},

  {{CA, CB, CG, OD1}, 1.231000, 120.800000, 0, 2},
  {{OD1, CB, CG, ND2}, 1.328000, 116.400000, 180.0, 0},

  {{CB, CG, ND2, HD21}, 1.020000, 120.000000, 0.0, 0},
  {{CB, CG, ND2, HD22}, 1.020000, 120.000000, 180.0, 0},
};

const cAminoResidue::residueRebuildingContext cAminoResidue::ASPresidueRebuildingData[] = {
  {{N, C, CA, HA}, 1.090000, 109.000000, -121.435662, 0},
  {{N, C, CA, CB}, 1.530000, 110.100000, +122.815261, 0},

  {{N, CA, CB, CG}, 1.516000, 112.600000, 0, 1},
  {{CG, CA, CB, HB1}, 1.090000, 109.000000, +119.823949, 0},
  {{CG, CA, CB, HB2}, 1.090000, 109.000000, -119.823949, 0},

  {{CA, CB, CG, OD1}, 1.249000, 118.400000, 0, 2},
  {{OD1, CB, CG, OD2}, 1.249000, 118.400000, 180.0, 0},
};

const cAminoResidue::residueRebuildingContext cAminoResidue::CYSresidueRebuildingData[] = {
  {{N, C, CA, HA}, 1.090000, 109.000000, -121.435662, 0},
  {{N, C, CA, CB}, 1.530000, 110.100000, +122.815261, 0},

  {{N, CA, CB, SG}, 1.808000, 114.400000, 0, 1},
  {{SG, CA, CB, HB1}, 1.090000, 109.000000, +121.002148, 0},
  {{SG, CA, CB, HB2}, 1.090000, 109.000000, -121.002148, 0},

  // {{CA, CB, SG, HG}, 1.300000, 109.000000, 0, 2},
  {{CA, CB, SG, HG}, 1.300000, 109.000000, 180.0, 0},
};

const cAminoResidue::residueRebuildingContext cAminoResidue::GLNresidueRebuildingData[] = {
  {{N, C, CA, HA}, 1.090000, 109.000000, -121.435662, 0},
  {{N, C, CA, CB}, 1.530000, 110.100000, +122.815261, 0},

  {{N, CA, CB, CG}, 1.520000, 114.100000, 0, 1},
  {{CG, CA, CB, HB1}, 1.090000, 109.000000, +120.800876, 0},
  {{CG, CA, CB, HB2}, 1.090000, 109.000000, -120.800876, 0},

  {{CA, CB, CG, CD}, 1.516000, 112.600000, 0, 2},
  {{CD, CB, CG, HG1}, 1.090000, 109.000000, +119.823949, 0},
  {{CD, CB, CG, HG2}, 1.090000, 109.000000, -119.823949, 0},

  {{CB, CG, CD, OE1}, 1.231000, 120.800000, 0, 3},
  {{OE1, CG, CD, NE2}, 1.328000, 116.400000, 180.0, 0},

  {{CG, CD, NE2, HE21}, 1.020000, 120.000000, 0.0, 0},
  {{CG, CD, NE2, HE22}, 1.020000, 120.000000, 180.0, 0},
};

const cAminoResidue::residueRebuildingContext cAminoResidue::GLUresidueRebuildingData[] = {
  {{N, C, CA, HA}, 1.090000, 109.000000, -121.435662, 0},
  {{N, C, CA, CB}, 1.530000, 110.100000, +122.815261, 0},

  {{N, CA, CB, CG}, 1.520000, 114.100000, 0, 1},
  {{CG, CA, CB, HB1}, 1.090000, 109.000000, +120.800876, 0},
  {{CG, CA, CB, HB2}, 1.090000, 109.000000, -120.800876, 0},

  {{CA, CB, CG, CD}, 1.516000, 112.600000, 0, 2},
  {{CD, CB, CG, HG1}, 1.090000, 109.000000, +119.823949, 0},
  {{CD, CB, CG, HG2}, 1.090000, 109.000000, -119.823949, 0},

  {{CB, CG, CD, OE1}, 1.249000, 118.400000, 0, 3},
  {{OE1, CG, CD, OE2}, 1.249000, 118.400000, 180.0, 0},
};

const cAminoResidue::residueRebuildingContext cAminoResidue::GLYresidueRebuildingData[] = {
  {{N, C, CA, HA1}, 1.090000, 109.000000, +122.286710, 0},
  {{N, C, CA, HA2}, 1.090000, 109.000000, -122.286710, 0},
};

const cAminoResidue::residueRebuildingContext cAminoResidue::HISresidueRebuildingData[] = {
  {{N, C, CA, HA}, 1.090000, 109.000000, -121.435662, 0},
  {{N, C, CA, CB}, 1.530000, 110.100000, +122.815261, 0},

  {{N, CA, CB, CG}, 1.497000, 113.800000, 0, 1},
  {{CG, CA, CB, HB1}, 1.090000, 109.000000, +120.601614, 0},
  {{CG, CA, CB, HB2}, 1.090000, 109.000000, -120.601614, 0},

  {{CA, CB, CG, ND1}, 1.378000, 122.700000, 0, 2},
  {{CB, CG, ND1, CE1}, 1.321000, 109.300000, 180.0, 0},
  {{CB, CG, ND1, HD1}, 1.020000, 125.350000, 0.0, 0},
  {{CG, ND1, CE1, NE2}, 1.321000, 108.400000, 0.0, 0},
  {{CG, ND1, CE1, HE1}, 1.080000, 125.800000, 180.0, 0},
  {{ND1, CE1, NE2, CD2}, 1.374000, 109.000000, 0.0, 0},
  {{ND1, CE1, NE2, HE2}, 1.020000, 125.500000, 180.0, 0},
  {{CE1, NE2, CD2, HD2}, 1.080000, 126.400000, 180.0, 0},
};

const cAminoResidue::residueRebuildingContext cAminoResidue::ILEresidueRebuildingData[] = {
  {{N, C, CA, HA}, 1.090000, 109.000000, -121.435662, 0},
  {{N, C, CA, CB}, 1.540000, 109.100000, +123.388988, 0},

  {{N, CA, CB, CG1}, 1.521000, 110.500000, 0, 1},
  {{CG1, CA, CB, CG2}, 1.530000, 110.400000, -122.797488, 0},
  {{CG1, CA, CB, HB}, 1.090000, 109.000000, +119.758560, 0},

  {{CA, CB, CG2, HG23}, 1.090000, 109.000000, 60.0, 0},
  {{CA, CB, CG2, HG22}, 1.090000, 109.000000, 180.0, 0},
  {{CA, CB, CG2, HG21}, 1.090000, 109.000000, 300.0, 0},

  {{CA, CB, CG1, CD1}, 1.513000, 113.800000, 0, 2},
  {{CD1, CB, CG1, HG11}, 1.090000, 109.000000, +120.601614, 0},
  {{CD1, CB, CG1, HG12}, 1.090000, 109.000000, -120.601614, 0},

  {{CB, CG1, CD1, HD13}, 1.090000, 109.000000, 60.0, 0},
  {{CB, CG1, CD1, HD12}, 1.090000, 109.000000, 180.0, 0},
  {{CB, CG1, CD1, HD11}, 1.090000, 109.000000, 300.0, 0},
};

const cAminoResidue::residueRebuildingContext cAminoResidue::LEUresidueRebuildingData[] = {
  {{N, C, CA, HA}, 1.090000, 109.000000, -121.435662, 0},
  {{N, C, CA, CB}, 1.530000, 110.100000, +122.815261, 0},

  {{N, CA, CB, CG}, 1.530000, 116.300000, 0, 1},
  {{CG, CA, CB, HB1}, 1.090000, 109.000000, +122.326073, 0},
  {{CG, CA, CB, HB2}, 1.090000, 109.000000, -122.326073, 0},

  {{CA, CB, CG, CD1}, 1.521000, 110.700000, 0, 2},
  {{CD1, CB, CG, CD2}, 1.521000, 110.700000, +123.270704, 0},
  {{CD1, CB, CG, HG}, 1.090000, 109.000000, -118.651949, 0},

  {{CB, CG, CD1, HD13}, 1.090000, 109.000000, 60, 0},
  {{CB, CG, CD1, HD12}, 1.090000, 109.000000, 180.0, 0},
  {{CB, CG, CD1, HD11}, 1.090000, 109.000000, 300.0, 0},

  {{CB, CG, CD2, HD23}, 1.090000, 109.000000, 60.0, 0},
  {{CB, CG, CD2, HD22}, 1.090000, 109.000000, 180.0, 0},
  {{CB, CG, CD2, HD21}, 1.090000, 109.000000, 300.0, 0},
};

const cAminoResidue::residueRebuildingContext cAminoResidue::LYSresidueRebuildingData[] = {
  {{N, C, CA, HA}, 1.090000, 109.000000, -121.435662, 0},
  {{N, C, CA, CB}, 1.530000, 110.100000, +122.815261, 0},

  {{N, CA, CB, CG}, 1.520000, 114.100000, 0, 1},
  {{CG, CA, CB, HB1}, 1.090000, 109.000000, +120.800876, 0},
  {{CG, CA, CB, HB2}, 1.090000, 109.000000, -120.800876, 0},

  {{CA, CB, CG, CD}, 1.520000, 111.300000, 0, 2},
  {{CD, CB, CG, HG1}, 1.090000, 109.000000, +119.014578, 0},
  {{CD, CB, CG, HG2}, 1.090000, 109.000000, -119.014578, 0},

  {{CB, CG, CD, CE}, 1.520000, 111.300000, 0, 3},
  {{CE, CG, CD, HD1}, 1.090000, 109.000000, +119.014578, 0},
  {{CE, CG, CD, HD2}, 1.090000, 109.000000, -119.014578, 0},

  {{CG, CD, CE, NZ}, 1.489000, 111.900000, 0, 4},
  {{NZ, CD, CE, HE1}, 1.090000, 109.000000, +119.384015, 0},
  {{NZ, CD, CE, HE2}, 1.090000, 109.000000, -119.384015, 0},

  {{CD, CE, NZ, HZ1}, 1.040000, 110.000000, 180.0, 0},
  {{CD, CE, NZ, HZ2}, 1.040000, 110.000000, 300.0, 0},
  {{CD, CE, NZ, HZ3}, 1.040000, 110.000000, 60.0, 0},
};

const cAminoResidue::residueRebuildingContext cAminoResidue::METresidueRebuildingData[] = {
  {{N, C, CA, HA}, 1.090000, 110.000000, -122.111265, 0},
  {{N, C, CA, CB}, 1.530000, 110.100000, +122.815261, 0},

  {{N, CA, CB, CG}, 1.520000, 114.100000, 0, 1},
  {{CG, CA, CB, HB1}, 1.090000, 109.000000, +120.800876, 0},
  {{CG, CA, CB, HB2}, 1.090000, 109.000000, -120.800876, 0},

  {{CA, CB, CG, SD}, 1.803000, 112.700000, 0, 2},
  {{SD, CB, CG, HG1}, 1.090000, 109.000000, +119.887603, 0},
  {{SD, CB, CG, HG2}, 1.090000, 109.000000, -119.887603, 0},

  {{CB, CG, SD, CE}, 1.791000, 100.900000, 0, 3},

  {{CG, SD, CE, HE1}, 1.090000, 109.000000, 180.0, 0},
  {{CG, SD, CE, HE2}, 1.090000, 109.000000, 300.0, 0},
  {{CG, SD, CE, HE3}, 1.090000, 109.000000, 60.0, 0},
};

//TODO: Check geometry
const cAminoResidue::residueRebuildingContext cAminoResidue::MSEresidueRebuildingData[] = {
  {{N, C, CA, HA}, 1.090000, 110.000000, -122.111265, 0},
  {{N, C, CA, CB}, 1.530000, 110.100000, +122.815261, 0},

  {{N, CA, CB, CG}, 1.520000, 114.100000, 0, 1},
  {{CG, CA, CB, HB1}, 1.090000, 109.000000, +120.800876, 0},
  {{CG, CA, CB, HB2}, 1.090000, 109.000000, -120.800876, 0},

  {{CA, CB, CG, SE}, 1.803000, 112.700000, 0, 2},
  {{SE, CB, CG, HG1}, 1.090000, 109.000000, +119.887603, 0},
  {{SE, CB, CG, HG2}, 1.090000, 109.000000, -119.887603, 0},

  {{CB, CG, SE, CE}, 1.791000, 100.900000, 0, 3},

  {{CG, SE, CE, HE1}, 1.090000, 109.000000, 180.0, 0},
  {{CG, SE, CE, HE2}, 1.090000, 109.000000, 300.0, 0},
  {{CG, SE, CE, HE3}, 1.090000, 109.000000, 60.0, 0},
};

const cAminoResidue::residueRebuildingContext cAminoResidue::PHEresidueRebuildingData[] = {
  {{N, C, CA, HA}, 1.090000, 109.000000, -121.435662, 0},
  {{N, C, CA, CB}, 1.530000, 110.100000, +122.815261, 0},

  {{N, CA, CB, CG}, 1.502000, 113.800000, 0, 1},
  {{CG, CA, CB, HB1}, 1.090000, 109.000000, +120.601614, 0},
  {{CG, CA, CB, HB2}, 1.090000, 109.000000, -120.601614, 0},

  {{CA, CB, CG, CD1}, 1.384000, 120.700000, 0, 2},
  {{CB, CG, CD1, CE1}, 1.382000, 120.700000, 180.0, 0},
  {{CB, CG, CD1, HD1}, 1.080000, 119.650000, 0.0, 0},
  {{CG, CD1, CE1, CZ}, 1.382000, 120.000000, 0.0, 0},
  {{CG, CD1, CE1, HE1}, 1.080000, 120.000000, 180.0, 0},
  {{CD1, CE1, CZ, CE2}, 1.382000, 120.000000, 0.0, 0},
  {{CD1, CE1, CZ, HZ}, 1.080000, 120.000000, 180.0, 0},
  {{CE1, CZ, CE2, CD2}, 1.382000, 120.000000, 0.0, 0},
  {{CE1, CZ, CE2, HE2}, 1.080000, 120.000000, 180.0, 0},
  {{CZ, CE2, CD2, HD2}, 1.080000, 119.650000, 180.0, 0},
};

const cAminoResidue::residueRebuildingContext cAminoResidue::PROresidueRebuildingData[] = {
  {{N, C, CA, HA}, 1.090000, 109.000000, -121.823876, 0},
  {{N, C, CA, CB}, 1.530000, 110.100000, +113.850977, 0},

  {{N, CA, CB, CG}, 1.492000, 104.500000, 0, 1},
  {{CG, CA, CB, HB1}, 1.090000, 109.000000, +117.560030, 0},
  {{CG, CA, CB, HB2}, 1.090000, 109.000000, -117.560030, 0},

  {{CA, CB, CG, CD}, 1.503000, 106.100000, 0, 2},
  {{CD, CB, CG, HG1}, 1.090000, 109.000000, +118.416600, 0},
  {{CD, CB, CG, HG2}, 1.090000, 109.000000, -118.416600, 0},

  // {{CB, CG, CD, N}, 1.473000, 103.200000, 0, 3},
  {{N, CG, CD, HD1}, 1.090000, 109.000000, +115.749241, 0},
  {{N, CG, CD, HD2}, 1.090000, 109.000000, -115.749241, 0},
};

const cAminoResidue::residueRebuildingContext cAminoResidue::SERresidueRebuildingData[] = {
  {{N, C, CA, HA}, 1.090000, 109.000000, -121.435662, 0},
  {{N, C, CA, CB}, 1.530000, 110.100000, +122.815261, 0},

  {{N, CA, CB, OG}, 1.417000, 111.100000, 0, 1},
  {{OG, CA, CB, HB1}, 1.090000, 109.000000, +120.128273, 0},
  {{OG, CA, CB, HB2}, 1.090000, 109.000000, -120.128273, 0},

  {{CA, CB, OG, HG}, 0.980000, 110.000000, 180.0, 0},
};

//TODO: Check geometry
const cAminoResidue::residueRebuildingContext cAminoResidue::SECresidueRebuildingData[] = {
  {{N, C, CA, HA}, 1.090000, 109.000000, -121.435662, 0},
  {{N, C, CA, CB}, 1.530000, 110.100000, +122.815261, 0},

  {{N, CA, CB, SE}, 1.808000, 114.400000, 0, 1},
  {{SE, CA, CB, HB1}, 1.090000, 109.000000, +121.002148, 0},
  {{SE, CA, CB, HB2}, 1.090000, 109.000000, -121.002148, 0},

  {{CA, CB, SE, HG}, 1.300000, 109.000000, 180.0, 0},
};

const cAminoResidue::residueRebuildingContext cAminoResidue::THRresidueRebuildingData[] = {
  {{N, C, CA, HA}, 1.090000, 109.000000, -121.435662, 0},
  {{N, C, CA, CB}, 1.540000, 109.100000, +123.388988, 0},

  {{N, CA, CB, OG1}, 1.521000, 110.500000, 0, 1},
  {{OG1, CA, CB, CG2}, 1.433000, 109.600000, -120.510600, 0},
  {{OG1, CA, CB, HB}, 1.090000, 109.000000, +118.532545, 0},

  {{CA, CB, OG1, HG1}, 0.980000, 110.000000, 180.0, 0},

  {{CA, CB, CG2, HG21}, 1.090000, 110.000000, 180.0, 0},
  {{CA, CB, CG2, HG22}, 1.090000, 110.000000, 300.0, 0},
  {{CA, CB, CG2, HG23}, 1.090000, 110.000000, 60.0, 0},
};

const cAminoResidue::residueRebuildingContext cAminoResidue::TRPresidueRebuildingData[] = {
  {{N, C, CA, HA}, 1.090000, 109.000000, -121.435662, 0},
  {{N, C, CA, CB}, 1.530000, 110.100000, +122.815261, 0},

  {{N, CA, CB, CG}, 1.498000, 113.600000, 0, 1},
  {{CG, CA, CB, HB1}, 1.090000, 109.000000, +120.469872, 0},
  {{CG, CA, CB, HB2}, 1.090000, 109.000000, -120.469872, 0},

  {{CA, CB, CG, CD1}, 1.365000, 126.900000, 0, 2},
  {{CB, CG, CD1, NE1}, 1.374000, 110.200000, 180.0, 0},
  {{CB, CG, CD1, HD1}, 1.080000, 124.900000, 0.0, 0},
  {{CG, CD1, NE1, CE2}, 1.370000, 108.900000, 0.0, 0},
  {{CG, CD1, NE1, HE1}, 1.020000, 125.550000, 180.0, 0},
  {{CD1, NE1, CE2, CD2}, 1.409000, 107.400000, 0.0, 0},
  {{NE1, CE2, CD2, CE3}, 1.398000, 118.800000, 180.0, 0},
  {{CE2, CD2, CE3, CZ3}, 1.382000, 118.600000, 0.0, 0},
  {{CE2, CD2, CE3, HE3}, 1.080000, 120.700000, 180.0, 0},
  {{CD2, CE3, CZ3, CH2}, 1.400000, 121.100000, 0.0, 0},
  {{CD2, CE3, CZ3, HZ3}, 1.080000, 119.450000, 180.0, 0},
  {{CE3, CZ3, CH2, CZ2}, 1.368000, 121.500000, 0.0, 0},
  {{CE3, CZ3, CH2, HH2}, 1.080000, 119.250000, 180.0, 0},
  {{CZ3, CH2, CZ2, CE2}, 1.394000, 117.500000, 0.0, 0},
  {{CZ3, CH2, CZ2, HZ2}, 1.080000, 121.250000, 180.0, 0},
};

const cAminoResidue::residueRebuildingContext cAminoResidue::TYRresidueRebuildingData[] = {
  {{N, C, CA, HA}, 1.090000, 109.000000, -121.435662, 0},
  {{N, C, CA, CB}, 1.530000, 110.100000, +122.815261, 0},

  {{N, CA, CB, CG}, 1.512000, 113.900000, 0, 1},
  {{CG, CA, CB, HB1}, 1.090000, 109.000000, +120.667814, 0},
  {{CG, CA, CB, HB2}, 1.090000, 109.000000, -120.667814, 0},

  {{CA, CB, CG, CD1}, 1.389000, 120.800000, 0, 2},
  {{CB, CG, CD1, CE1}, 1.382000, 121.200000, 180.0, 0},
  {{CB, CG, CD1, HD1}, 1.080000, 119.400000, 0.0, 0},
  {{CG, CD1, CE1, CZ}, 1.378000, 119.600000, 0.0, 0},
  {{CG, CD1, CE1, HE1}, 1.080000, 120.200000, 180.0, 0},
  {{CD1, CE1, CZ, CE2}, 1.378000, 120.300000, 0.0, 0},
  {{CD1, CE1, CZ, OH}, 1.376000, 119.900000, 180.0, 0},
  {{CE1, CZ, OH, HH}, 0.980000, 110.000000, 180.0, 0},
  {{CE1, CZ, CE2, CD2}, 1.382000, 119.600000, 0.0, 0},
  {{CE1, CZ, CE2, HE2}, 1.080000, 120.200000, 180.0, 0},
  {{CZ, CE2, CD2, HD2}, 1.080000, 119.400000, 180.0, 0},
};

const cAminoResidue::residueRebuildingContext cAminoResidue::VALresidueRebuildingData[] = {
  {{N, C, CA, HA}, 1.090000, 109.000000, -121.435662, 0},
  {{N, C, CA, CB}, 1.540000, 109.100000, +123.388988, 0},

  {{N, CA, CB, CG1}, 1.521000, 110.400000, 0, 1},
  {{CG1, CA, CB, CG2}, 1.521000, 110.400000, +122.855892, 0},
  {{CG1, CA, CB, HB}, 1.090000, 109.000000, -118.473114, 0},

  {{CA, CB, CG1, HG11}, 1.090000, 109.000000, 180.0, 0},
  {{CA, CB, CG1, HG11}, 1.090000, 109.000000, 300.0, 0},
  {{CA, CB, CG1, HG11}, 1.090000, 109.000000, 60.0, 0},

  {{CA, CB, CG2, HG21}, 1.090000, 109.000000, 180.0, 0},
  {{CA, CB, CG2, HG22}, 1.090000, 109.000000, 300.0, 0},
  {{CA, CB, CG2, HG23}, 1.090000, 109.000000, 60.0, 0},
};

cAminoResidue::cAminoResidue(const std::string &resName_, char chainId_,
                                   const std::string &segId_, size_t seqNumber_, cResidue* previous_)
    : cResidue(resName_, chainId_, segId_, seqNumber_, previous_) {

  std::map<std::string, eAminoType>::const_iterator it = mAminoType.find(resName);
  assert(it != cAminoResidue::mAminoType.end());
  aminoType_ = it->second;

#define DEFINE_ATOMS(...) { static eAtomLabel cl[] = {__VA_ARGS__}; \
                            for (size_t i = 0; i < sizeof(cl) / sizeof(eAtomLabel); ++i) \
                              atom_map[cl[i]] = NULL; }
#define DEFINE_REBUILDING_DATA(__data__) { \
      rotamersContext.residueRebuildingDataSize = sizeof(__data__) / sizeof(residueRebuildingContext); \
      rotamersContext.residueRebuildingData = new residueRebuildingContext[rotamersContext.residueRebuildingDataSize]; \
      memcpy(rotamersContext.residueRebuildingData, __data__, sizeof(__data__)); }
#define DEFINE_ALIASE(al_1, al_2) mAtomAliases[al_1] = al_2;

  switch (aminoType_) {
    case ALA:
      DEFINE_ATOMS( C, H, N, O, CA, HA, CB, HB3, HB2, HB1 );
      DEFINE_REBUILDING_DATA(ALAresidueRebuildingData);
      DEFINE_ALIASE(HA2, HA);
      break;
    case ARG:
      DEFINE_ATOMS( C, H, N, O, CA, HA, CB, HB2, HB1, CD, HD2, HD1, HE,
                    NE, CG, HG2, HG1, NH1, NH2, HH11, HH21, HH12, HH22, CZ );
      DEFINE_REBUILDING_DATA(ARGresidueRebuildingData);
      DEFINE_ALIASE(HA2, HA);
      DEFINE_ALIASE(HB3, HB1);
      DEFINE_ALIASE(HG3, HG1);
      DEFINE_ALIASE(HD3, HD1);
      break;
    case ASN:
      DEFINE_ATOMS( C, N, H, O, CA, CB, ND2, OD1, CG, HA, HB2, HB1, HD22, HD21 );
      DEFINE_REBUILDING_DATA(ASNresidueRebuildingData);
      DEFINE_ALIASE(HA2, HA);
      DEFINE_ALIASE(HB3, HB1);
      break;
    case ASP:
      DEFINE_ATOMS( C, N, H, O, CA, CB, OD1, OD2, CG, HA, HB2, HB1 );
      DEFINE_REBUILDING_DATA(ASPresidueRebuildingData);
      DEFINE_ALIASE(HA2, HA);
      DEFINE_ALIASE(HB3, HB1);
      break;
    case CYS:
      DEFINE_ATOMS( C, H, N, O, CA, HA, CB, HB2, HB1, HG, SG );
      DEFINE_REBUILDING_DATA(CYSresidueRebuildingData);
      DEFINE_ALIASE(HA2, HA);
      DEFINE_ALIASE(HG1, HG);
      DEFINE_ALIASE(HB3, HB1);
      break;
    case GLN:
      DEFINE_ATOMS( C, N, H, O, CA, CB, CD, NE2, OE1, CG, HA, HB2, HB1, HG2, HG1, HE22, HE21 );
      DEFINE_REBUILDING_DATA(GLNresidueRebuildingData);
      DEFINE_ALIASE(HA2, HA);
      DEFINE_ALIASE(HB3, HB1);
      DEFINE_ALIASE(HG3, HG1);
      break;
    case GLU:
      DEFINE_ATOMS( C, N, H, O, CA, CB, CD, OE1, OE2, CG, HA, HB2, HB1, HG2, HG1 );
      DEFINE_REBUILDING_DATA(GLUresidueRebuildingData);
      DEFINE_ALIASE(HA2, HA);
      DEFINE_ALIASE(HB3, HB1);
      DEFINE_ALIASE(HG3, HG1);
      break;
    case GLY:
      DEFINE_ATOMS( N, CA, C, O, H, HA2, HA1 );
      DEFINE_REBUILDING_DATA(GLYresidueRebuildingData);
      DEFINE_ALIASE(HA3, HA1);
      DEFINE_ALIASE(HA, HA1);
      break;
    case HIS:
      DEFINE_ATOMS( C, H, N, O, CA, HA, CB, HB2, HB1, CD2, HD1, HD2, ND1, CE1, HE1, HE2, NE2, CG );
      DEFINE_REBUILDING_DATA(HISresidueRebuildingData);
      DEFINE_ALIASE(HA2, HA);
      DEFINE_ALIASE(HB3, HB1);
      DEFINE_ALIASE(HD3, HD1);
      DEFINE_ALIASE(HE3, HE1);
      break;
    case ILE:
      DEFINE_ATOMS( C, H, N, O, CA, HA, CB, HB, CD1, HD13, HD12,
                    HD11, CG1, CG2, HG12, HG23, HG11, HG22, HG21 );
      DEFINE_REBUILDING_DATA(ILEresidueRebuildingData);
      DEFINE_ALIASE(HA2, HA);
      DEFINE_ALIASE(HG13, HG11);
      DEFINE_ALIASE(CD, CD1);
      DEFINE_ALIASE(HD1, HD11);
      DEFINE_ALIASE(HD2, HD12);
      DEFINE_ALIASE(HD3, HD13);
      break;
    case LEU:
      DEFINE_ATOMS( C, H, N, O, CA, HA, CB, HB2, HB1, CD1, CD2,
                    HD13, HD23, HD12, HD22, HD11, HD21, CG, HG );
      DEFINE_REBUILDING_DATA(LEUresidueRebuildingData);
      DEFINE_ALIASE(HA2, HA);
      DEFINE_ALIASE(HB3, HB1);
      break;
    case LYS:
      DEFINE_ATOMS( C, H, N, O, CA, HA, CB, HB2, HB1, CD, HD2, HD1,
                    CE, HE2, HE1, CG, HG2, HG1, NZ, HZ3, HZ2, HZ1 );
      DEFINE_REBUILDING_DATA(LYSresidueRebuildingData);
      DEFINE_ALIASE(HA2, HA);
      DEFINE_ALIASE(HB3, HB1);
      DEFINE_ALIASE(HG3, HG1);
      DEFINE_ALIASE(HD3, HD1);
      DEFINE_ALIASE(HE3, HE1);
      break;
    case MET:
      DEFINE_ATOMS( C, H, N, O, CA, HA, CB, HB2, HB1, SD, CE, HE3, HE2, HE1, CG, HG2, HG1 );
      DEFINE_REBUILDING_DATA(METresidueRebuildingData);
      DEFINE_ALIASE(HA2, HA);
      DEFINE_ALIASE(HG3, HG1);
      DEFINE_ALIASE(HB3, HB1);
      break;
    case MSE:
      DEFINE_ATOMS( C, H, N, O, CA, HA, CB, HB2, HB1, SE, CE, HE3, HE2, HE1, CG, HG2, HG1 );
      DEFINE_REBUILDING_DATA(MSEresidueRebuildingData);
      DEFINE_ALIASE(HA2, HA);
      DEFINE_ALIASE(HG3, HG1);
      DEFINE_ALIASE(HB3, HB1);
      break;
    case PHE:
      DEFINE_ATOMS( C, N, H, O, CA, CB, CD1, CD2, CE1, CE2, CG,
                    CZ, HA, HB2, HB1, HD1, HD2, HE1, HE2, HZ );
      DEFINE_REBUILDING_DATA(PHEresidueRebuildingData);
      DEFINE_ALIASE(HA2, HA);
      DEFINE_ALIASE(HB3, HB1);
      break;
    case PRO:
      DEFINE_ATOMS( C, CA, CB, CD, CG, N, O, HA, HB2, HB1, HG2, HG1, HD2, HD1 );
      DEFINE_REBUILDING_DATA(PROresidueRebuildingData);
      DEFINE_ALIASE(HA2, HA);
      DEFINE_ALIASE(HB3, HB1);
      DEFINE_ALIASE(HD3, HD1);
      DEFINE_ALIASE(HG3, HG1);
      break;
    case SER:
      DEFINE_ATOMS( C, H, N, O, CA, HA, CB, HB2, HB1, HG, OG );
      DEFINE_REBUILDING_DATA(SERresidueRebuildingData);
      DEFINE_ALIASE(HA2, HA);
      DEFINE_ALIASE(HG1, HG);
      DEFINE_ALIASE(HB3, HB1);
      break;
    case SEC:
      DEFINE_ATOMS( C, H, N, O, CA, HA, CB, HB2, HB1, HG, SE );
      DEFINE_REBUILDING_DATA(SECresidueRebuildingData);
      DEFINE_ALIASE(HA2, HA);
      DEFINE_ALIASE(HG1, HG);
      DEFINE_ALIASE(HB3, HB1);
      break;
    case THR:
      DEFINE_ATOMS( C, H, N, O, CA, HA, CB, HB, CG2, HG1, OG1, HG23, HG22, HG21 );
      DEFINE_REBUILDING_DATA(THRresidueRebuildingData);
      DEFINE_ALIASE(HA2, HA);
      break;
    case TRP:
      DEFINE_ATOMS( C, H, N, O, CA, HA, CB, HB2, HB1, CD1, CD2, HD1, CE2,
                    CE3, HE1, HE3, NE1, CG, CH2, HH2, CZ2, CZ3, HZ2, HZ3 );
      DEFINE_REBUILDING_DATA(TRPresidueRebuildingData);
      DEFINE_ALIASE(HA2, HA);
      DEFINE_ALIASE(HB3, HB1);
      DEFINE_ALIASE(HZ1, HZ3);
      break;
    case TYR:
      DEFINE_ATOMS( C, N, H, O, CA, CB, CD1, CD2, CE1, CE2, CG,
                    OH, CZ, HA, HB2, HB1, HD1, HD2, HE1, HE2, HH );
      DEFINE_REBUILDING_DATA(TYRresidueRebuildingData);
      DEFINE_ALIASE(HA2, HA);
      DEFINE_ALIASE(HB3, HB1);
      DEFINE_ALIASE(HD3, HD1);
      DEFINE_ALIASE(HE3, HE1);
      break;
    case VAL:
      DEFINE_ATOMS( C, H, N, O, CA, HA, CB, HB, CG1, CG2, HG13, HG23, HG12, HG22, HG11, HG21 );
      DEFINE_REBUILDING_DATA(VALresidueRebuildingData);
      DEFINE_ALIASE(HA2, HA);
      break;
    default:
      assert(false);
  }
  DEFINE_ATOMS( HT1, HT2, HT3, OXT );
  DEFINE_ALIASE(O1, O);
  DEFINE_ALIASE(O2, OXT);
}

cAminoResidue::eAtomLabel cAminoResidue::get_atom_label(const cAtom &atom) const {
  std::map<std::string, eAtomLabel>::const_iterator it = mAtomLabel.find(atom.name);
  if (it == mAtomLabel.end()) {
    return ATOM_LABEL_END;
  }
  eAtomLabel al = it->second;
  auto alias_it = mAtomAliases.find(al);
  if (alias_it != mAtomAliases.end()) {
    al = alias_it->second;
  }
  return al;
}

void cAminoResidue::addAtom(cAtom *atom) {
  assert(atom);

  cResidue::addAtom(atom);
  assert(numAtoms());
  if (atom != atoms_.back())
    return;

  auto al = get_atom_label(*atom);

  if (al == ATOM_LABEL_END) {
    cResidue::removeAtom(std::prev(atoms_.end()));
    throw std::runtime_error(std::string("Unknown atom ") + atom->name);
  }

  std::map<eAtomLabel, cAtom *>::iterator it_ = atom_map.find(al);
  if (it_ == atom_map.end()) {
    cResidue::removeAtom(std::prev(atoms_.end()));
    throw std::runtime_error(std::string("Amino acid ") + resName +
                             " doesn't have atom " + atom->name);
  }
  if (it_->second) {
    std::ostringstream string_stream;
    string_stream << "Redefinition of atom " << atom->name << "\n"
                  << *it_->second
                  << *atom;
    std::string error_info = string_stream.str();
    if (error_info.back() == '\n')
      error_info.pop_back();

    cResidue::removeAtom(std::prev(atoms_.end()));
    throw std::runtime_error(error_info);
  }
  it_->second = atom;

  // CHARMM Params
  cCharmmParams::setCharmmParams(atom, it_->first);
}

void cAminoResidue::removeAtom(std::vector<cAtom *>::iterator atom_it) {
  auto atom = *atom_it;
  auto it = mAtomLabel.find(atom->name);
  assert(it != mAtomLabel.end());

  eAtomLabel al = it->second;
  auto alias_it = mAtomAliases.find(al);
  if (alias_it != mAtomAliases.end())
    al = alias_it->second;

  auto it_ = atom_map.find(al);
  assert(it_ != atom_map.end());
  assert(it_->second);
  it_->second = nullptr;

  cResidue::removeAtom(atom_it);
}

void cAminoResidue::predictSideChain(const cRotamerLibrary &confLib) {
  const std::vector<cRotamer> *pConf = confLib.getRotamers(*this);
  if (pConf)
    rotamersContext.rotamers = *pConf;
}


bool cAminoResidue::isBackboneValid() const {
  const cAtom *CA = atom_map.find(eAtomLabel::CA)->second,
                  *N = atom_map.find(eAtomLabel::N)->second,
                  *C = atom_map.find(eAtomLabel::C)->second;
  if (!N || !CA || !C)
    return false;
  return !((C->getPosition() - CA->getPosition())
          ^ (N->getPosition() - CA->getPosition())).isZero();
}


void setDihedral(const cAtom *A, const cAtom *B, const cAtom *C,
                 cAtom *D, double CD, double B_C_D, double A_B_C_D) {

  if (!A || !B || !C || !D)
    return;

  cVector3 D_pos;
  math_functions::setDihedral(A->getPosition(), B->getPosition(), C->getPosition(),
                              &D_pos, CD, B_C_D, A_B_C_D);
  D->setPosition(D_pos);
}

bool getDist(const cAtom *A, const cAtom *B, double *dist) {
  if (!A || !B)
    return false;
  *dist = (A->getPosition() - B->getPosition()).norm();
  return true;
}

bool getAngle(const cAtom *A, const cAtom *B, const cAtom *C, double *angle) {
  if (!A || !B || !C)
    return false;

  return math_functions::computeAngle(A->getPosition(),
                                      B->getPosition(),
                                      C->getPosition(), angle);
}

bool getDihedral(const cAtom *A, const cAtom *B,
                 const cAtom *C, const cAtom *D, double *dihedral) {
  if (!A || !B || !C || !D)
    return false;

  return math_functions::computeDihedralAngle(A->getPosition(), B->getPosition(),
                                              C->getPosition(), D->getPosition(),
                                              dihedral);
}

void cAminoResidue::handlePhiPsi(cAminoResidue *nextRes) {
  cAtom *N0 = atom_map[N],
           *CA0 = atom_map[CA],
           *C0 = atom_map[C],
           *N1 = nextRes->atom_map[N],
           *CA1 = nextRes->atom_map[CA],
           *C1 = nextRes->atom_map[C];

  if (!N0 || !C0 || !N1 || !C1)
    return;

  double C0_N1 = (C0->getPosition() - N1->getPosition()).norm();
  double N0_C1 = (C1->getPosition() - N0->getPosition()).norm();
  if (C0_N1 < N0_C1) {
    if (C0_N1 < 2.0) {
      getDihedral(N0, CA0, C0, N1, &psi);
      getDihedral(C0, N1, CA1, C1, &nextRes->phi);
    }
  } else {
    if (N0_C1 < 2.0) {
      getDihedral(N1, CA1, C1, N0, &nextRes->psi);
      getDihedral(C1, N0, CA0, C0, &phi);
    }
  }
}

void cAminoResidue::addAtom(eAtomLabel atomLabel, bool withHydrogens) {
  std::map<std::string, eAtomLabel>::const_iterator it;
  for (it = mAtomLabel.begin(); it != mAtomLabel.end(); ++it) {
    if (it->second == atomLabel && it->first[0] >= 'A')
      break;
  }
  assert(it != mAtomLabel.end());

  if (!withHydrogens && it->first[0] == 'H')
    return;

  char buffer[81];

  sprintf(buffer, "ATOM  %5d %-4s %3s %c%4zd    %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s %c  ",
                  0, it->first.c_str(), resName.c_str(), chainId, seqNumber,
                  0.0, 0.0, 0.0, 1.0, 0.0, segId.c_str(), it->first[0]);

  cAtom *atom = new cAtom();
  if (!atom->parseFromPDBField(buffer)) {
    assert(false && "Can't parse atom PDB line");
  }
  addAtom(atom);
}

void cAminoResidue::complete(bool withHydrogens) {
  for (size_t i = 0; i < rotamersContext.residueRebuildingDataSize; ++i) {
    cAtom *A = atom_map[rotamersContext.residueRebuildingData[i].A[0]],
             *B = atom_map[rotamersContext.residueRebuildingData[i].A[1]],
             *C = atom_map[rotamersContext.residueRebuildingData[i].A[2]],
             *D = atom_map[rotamersContext.residueRebuildingData[i].A[3]];
    if (!A || !B || !C || D)
      continue;

    double dist = rotamersContext.residueRebuildingData[i].dist,
           angle = rotamersContext.residueRebuildingData[i].angle,
           dihedral = rotamersContext.residueRebuildingData[i].dihedral;

    if (rotamersContext.residueRebuildingData[i].chi > 0) {
      if (!numRotamers())
        continue;

      dihedral = rotamersContext.rotamers[0].chi[rotamersContext.residueRebuildingData[i].chi - 1];
    }

    addAtom(rotamersContext.residueRebuildingData[i].A[3], withHydrogens);

    D = atom_map[rotamersContext.residueRebuildingData[i].A[3]];
    // Compute coordinates of the atom D
    setDihedral(A, B, C, D, dist, angle, dihedral);
  }
}

void cAminoResidue::fitResidueGeometry() {
  for (size_t i = 0; i < rotamersContext.residueRebuildingDataSize; ++i) {
    cAtom *A = atom_map[rotamersContext.residueRebuildingData[i].A[0]],
             *B = atom_map[rotamersContext.residueRebuildingData[i].A[1]],
             *C = atom_map[rotamersContext.residueRebuildingData[i].A[2]],
             *D = atom_map[rotamersContext.residueRebuildingData[i].A[3]];
    getDist(C, D, &rotamersContext.residueRebuildingData[i].dist);
    getAngle(B, C, D, &rotamersContext.residueRebuildingData[i].angle);
    getDihedral(A, B, C, D, &rotamersContext.residueRebuildingData[i].dihedral);
  }
}

#ifdef RASP_ENERGY
//Zhichao Miao, Yang Cao and Taijiao Jiang
const std::map<cAminoResidue::eAminoType, double> mScalingEnergyFactor = {
  {cAminoResidue::eAminoType::ALA, 0.00},
  {cAminoResidue::eAminoType::ARG, 1.40},
  {cAminoResidue::eAminoType::ASN, 0.90},
  {cAminoResidue::eAminoType::ASP, 1.30},
  {cAminoResidue::eAminoType::CYS, 1.00},
  {cAminoResidue::eAminoType::GLN, 1.10},
  {cAminoResidue::eAminoType::GLU, 1.20},
  {cAminoResidue::eAminoType::GLY, 0.00},
  {cAminoResidue::eAminoType::HIS, 1.90},
  {cAminoResidue::eAminoType::ILE, 1.00},
  {cAminoResidue::eAminoType::LEU, 1.10},
  {cAminoResidue::eAminoType::LYS, 1.00},
  {cAminoResidue::eAminoType::MET, 1.10},
  {cAminoResidue::eAminoType::PHE, 2.10},
  {cAminoResidue::eAminoType::PRO, 0.50},
  {cAminoResidue::eAminoType::SER, 0.90},
  {cAminoResidue::eAminoType::THR, 1.50},
  {cAminoResidue::eAminoType::TRP, 1.00},
  {cAminoResidue::eAminoType::TYR, 1.00},
  {cAminoResidue::eAminoType::VAL, 1.00},
};
#elif SCWRL4_ENERGY
//Georgii G. Krivov, Maxim V. Shapovalov and Roland L. Dunbrack, Jr
const std::map<cAminoResidue::eAminoType, double> mScalingEnergyFactor = {
  {cAminoResidue::eAminoType::ALA, 0.00},
  {cAminoResidue::eAminoType::ARG, 2.27},
  {cAminoResidue::eAminoType::ASN, 1.80},
  {cAminoResidue::eAminoType::ASP, 2.44},
  {cAminoResidue::eAminoType::CYS, 4.07},
  {cAminoResidue::eAminoType::GLN, 1.61},
  {cAminoResidue::eAminoType::GLU, 1.85},
  {cAminoResidue::eAminoType::GLY, 0.00},
  {cAminoResidue::eAminoType::HIS, 2.01},
  {cAminoResidue::eAminoType::ILE, 2.18},
  {cAminoResidue::eAminoType::LEU, 2.25},
  {cAminoResidue::eAminoType::LYS, 2.13},
  {cAminoResidue::eAminoType::MET, 1.95},
  {cAminoResidue::eAminoType::PHE, 1.71},
  {cAminoResidue::eAminoType::PRO, 0.76},
  {cAminoResidue::eAminoType::SER, 2.78},
  {cAminoResidue::eAminoType::THR, 2.96},
  {cAminoResidue::eAminoType::TRP, 3.24},
  {cAminoResidue::eAminoType::TYR, 2.00},
  {cAminoResidue::eAminoType::VAL, 1.62},
};
#endif

double cAminoResidue::setRotamer(size_t j, double *libEnergy) {
  assert(j < numRotamers());

  for (size_t i = 0; i < rotamersContext.residueRebuildingDataSize; ++i) {
    double dist = rotamersContext.residueRebuildingData[i].dist,
           angle = rotamersContext.residueRebuildingData[i].angle,
           dihedral = rotamersContext.residueRebuildingData[i].dihedral;
    cAtom *A = atom_map[rotamersContext.residueRebuildingData[i].A[0]],
             *B = atom_map[rotamersContext.residueRebuildingData[i].A[1]],
             *C = atom_map[rotamersContext.residueRebuildingData[i].A[2]],
             *D = atom_map[rotamersContext.residueRebuildingData[i].A[3]];

    if (rotamersContext.residueRebuildingData[i].chi > 0)
      dihedral = rotamersContext.rotamers[j].chi[rotamersContext.residueRebuildingData[i].chi - 1];

    setDihedral(A, B, C, D, dist, angle, dihedral);
  }
  if (libEnergy) {
    auto it = mScalingEnergyFactor.find(aminoType_);
    *libEnergy = std::min(51.0, -it->second * log(rotamersContext.rotamers[j].prob
                                                  / rotamersContext.rotamers[0].prob));
  }
  return rotamersContext.rotamers[j].prob;
}

double cAminoResidue::getChi(size_t n) const {
  assert(n);
  for (size_t i = 0; i < rotamersContext.residueRebuildingDataSize; ++i) {
    if (rotamersContext.residueRebuildingData[i].chi != n)
      continue;

    double chi_n;
    std::map<eAtomLabel, cAtom *>::const_iterator
      A = atom_map.find(rotamersContext.residueRebuildingData[i].A[0]),
      B = atom_map.find(rotamersContext.residueRebuildingData[i].A[1]),
      C = atom_map.find(rotamersContext.residueRebuildingData[i].A[2]),
      D = atom_map.find(rotamersContext.residueRebuildingData[i].A[3]);
    assert(A != atom_map.end());
    assert(B != atom_map.end());
    assert(C != atom_map.end());
    assert(D != atom_map.end());
    if (!A->second || !B->second || !C->second || !D->second)
      return 0.0;
    if (!getDihedral(A->second, B->second, C->second, D->second, &chi_n)) {
      assert(false && "Can't get chi angle");
    }
    return chi_n;
  }
  return 0.0;
  fprintf(stderr, "ERROR!!! Amino acid %s doesn't have angle chi_%zd\n", resName.c_str(), n);
  exit(1);
}

size_t cAminoResidue::numRotamers() const {
  return rotamersContext.rotamers.size();
}

const cAtom* cAminoResidue::getAtomByLabel(eAtomLabel atomLabel) const {
  std::map<eAtomLabel, cAtom *>::const_iterator it = atom_map.find(atomLabel);
  assert(it != atom_map.end() && "Unknown atom for the amino residue");
  // if (!it->second) {
  //   fprintf(stderr, "ERROR!!! Access to undefined atom in residue:\n%3s %c%4zd      %-4s\n",
  //                   resName.c_str(), chainId, seqNumber, segId.c_str());
  // }
  return it->second;
}

const std::map<std::string, cAminoResidue::eAtomLabel> cAminoResidue::mAtomLabel = {
  {"C", C},
  {"N", N},
  {"O", O},
  {"H", H}, {"HN", H},
  {"H1", HT1}, {"1H", HT1}, {"HT1", HT1}, {"1HT", HT1},
  {"H2", HT2}, {"2H", HT2}, {"HT2", HT2}, {"2HT", HT2},
  {"H3", HT3}, {"3H", HT3}, {"HT3", HT3}, {"3HT", HT3},

  {"CA", CA},
  {"HA", HA},
  {"HA1", HA1}, {"1HA", HA1},
  {"HA2", HA2}, {"2HA", HA2},
  {"HA3", HA3}, {"3HA", HA3},

  {"CB", CB},
  {"HB", HB},
  {"HB1", HB1}, {"1HB", HB1},
  {"HB2", HB2}, {"2HB", HB2},
  {"HB3", HB3}, {"3HB", HB3},

  {"CG", CG},
  {"CG1", CG1}, {"1CG", CG1},
  {"CG2", CG2}, {"2CG", CG2},
  {"HG", HG},
  {"HG1", HG1}, {"1HG", HG1},
  {"HG2", HG2}, {"2HG", HG2},
  {"HG3", HG3}, {"3HG", HG3},

  {"HG11", HG11}, {"1HG1", HG11},
  {"HG12", HG12}, {"2HG1", HG12},
  {"HG13", HG13}, {"3HG1", HG13},

  {"HG21", HG21}, {"1HG2", HG21},
  {"HG22", HG22}, {"2HG2", HG22},
  {"HG23", HG23}, {"3HG2", HG23},

  {"OG", OG},
  {"OG1", OG1}, {"1OG", OG1},
  {"SG", SG},

  {"CD", CD},
  {"CD1", CD1}, {"1CD", CD1},
  {"CD2", CD2}, {"2CD", CD2},
  {"HD", HD},
  {"HD1", HD1}, {"1HD", HD1},
  {"HD2", HD2}, {"2HD", HD2},
  {"HD3", HD3}, {"3HD", HD3},

  {"HD11", HD11}, {"1HD1", HD11},
  {"HD12", HD12}, {"2HD1", HD12},
  {"HD13", HD13}, {"3HD1", HD13},
  {"HD21", HD21}, {"1HD2", HD21},
  {"HD22", HD22}, {"2HD2", HD22},
  {"HD23", HD23}, {"3HD2", HD23},

  {"SD", SD},
  {"SE", SE},
  {"OD1", OD1}, {"1OD", OD1},
  {"OD2", OD2}, {"2OD", OD2},
  {"ND1", ND1}, {"1ND", ND1},
  {"ND2", ND2}, {"2ND", ND2},

  {"CE", CE},
  {"CE1", CE1}, {"1CE", CE1},
  {"CE2", CE2}, {"2CE", CE2},
  {"CE3", CE3}, {"3CE", CE3},
  {"HE", HE},
  {"HE1", HE1}, {"1HE", HE1},
  {"HE2", HE2}, {"2HE", HE2},
  {"HE3", HE3}, {"3HE", HE3},

  {"HE21", HE21}, {"1HE2", HE21},
  {"HE22", HE22}, {"2HE2", HE22},

  {"OE1", OE1}, {"1OE", OE1},
  {"OE2", OE2}, {"2OE", OE2},
  {"NE", NE},
  {"NE1", NE1}, {"1NE", NE1},
  {"NE2", NE2}, {"2NE", NE2},

  {"CZ", CZ},
  {"CZ2", CZ2}, {"2CZ", CZ2},
  {"CZ3", CZ3}, {"3CZ", CZ3},
  {"HZ", HZ},
  {"HZ1", HZ1}, {"1HZ", HZ1},
  {"HZ2", HZ2}, {"2HZ", HZ2},
  {"HZ3", HZ3}, {"3HZ", HZ3},
  {"NZ", NZ},

  {"CH2", CH2}, {"2CH", CH2},
  {"CH3", CH3}, {"3CH", CH3},
  {"NH1", NH1}, {"1NH", NH1},
  {"NH2", NH2}, {"2NH", NH2},
  {"OH", OH},
  {"HH", HH},
  {"HH11", HH11}, {"1HH1", HH11},

  {"HH12", HH12}, {"2HH1", HH12},
  {"HH2", HH2}, {"2HH", HH2},
  {"HH21", HH21}, {"1HH2", HH21},
  {"HH22", HH22}, {"2HH2", HH22},

  {"HH31", HH31}, {"1HH3", HH31},
  {"HH32", HH32}, {"2HH3", HH32},
  {"HH33", HH33}, {"3HH3", HH33},
  {"O1", O1}, {"O2", O2},
  {"OC1", O1}, {"OC2", O2},
  {"OT1", O}, {"OT2", OXT},
  {"OXT", OXT}

};

const std::map<std::string, cAminoResidue::eAminoType> cAminoResidue::mAminoType = {
  {"ALA", ALA},
  {"ARG", ARG},
  {"ASN", ASN},
  {"ASP", ASP},
  {"CYS", CYS},
  {"GLN", GLN},
  {"GLU", GLU},
  {"GLY", GLY},
  {"HIS", HIS},
  {"ILE", ILE},
  {"LEU", LEU},
  {"LYS", LYS},
  {"MET", MET},
  {"PHE", PHE},
  {"PRO", PRO},
  {"SER", SER},
  {"THR", THR},
  {"TRP", TRP},
  {"TYR", TYR},
  {"VAL", VAL},
  {"SEC", SEC},
  {"MSE", MSE},
};

//this is the version for oriented graphs
bool cAminoResidue::setLocalFrame(bool orient) {
    if (!this->getAtom("N")) return false;
    if (!this->getAtom("CA")) return false;
    if (!this->getAtom("C")) return false;
    origin = this->getAtom("CA")->getPosition();
    if (!orient) {
        vX = cVector3(1, 0, 0);
        vY = cVector3(0, 1, 0);
        vZ = cVector3(0, 0, 1);

        localToWorld = cSpatialTransform( cMatrix33(vX, vY, vZ), origin);
        worldToLocal = localToWorld.getInverse();

        return true;
    }

    // vector X is from CA to N
    vX = this->getAtom("N")->getPosition() - origin;
    vX.normalize();

    // vector Y is from CA to C
    vY = this->getAtom("C")->getPosition() - origin;
    vY -= vX * (vY | vX);
    vY.normalize();

    //vZ is orthogonal to the XY plane
    vZ = vX^vY;

    localToWorld = cSpatialTransform( cMatrix33(vX, vY, vZ), origin);
    worldToLocal = localToWorld.getInverse();
    return true;
}
