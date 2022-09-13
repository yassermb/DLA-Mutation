#pragma once

#include "cProtein.hpp"


namespace energy {

// Energy of the interaction of the side-chain of |res|
// with the neighbouring backbone atoms
double frameEnergy(const cResidue &res);

// Energy of the interaction of the side-chains of |res_1| and |res_2|
double pairEnergy(const cResidue &res_1, const cResidue &res_2);

// Total energy of the interaction of side-chain
// of residue |res| with the neighbouring atoms
double totalEnergy(const cResidue &res);

// RASP method of detecting |residues| -- amino residues that have effective
// contact among all amino residues in protein, that have more than 1 retamers.
// Neighbours for sidechain atoms from |residues| will be initialized among
// all other atoms (backbone) and for backbone atoms will be cleared
// return vector: neighbours for amino residues |residues|
// |residues| -- output parameter
std::vector<std::vector<cAminoResidue *>>
    initializeNeighbourResidues(cProtein *protein, std::vector<cAminoResidue *> *residues);

// Initialize neighbours for all atoms in protein
void initializeNeighbourAtoms(cProtein *protein);

// Number of atom pairs in clash
size_t numClashes(const cProtein &protein);

} // energy
