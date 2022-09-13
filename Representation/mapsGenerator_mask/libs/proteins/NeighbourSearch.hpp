#pragma once

#include <map>
#include <set>
#include <vector>
#include <algorithm>
#include <cstring>

#include "cAminoResidue.hpp"
#include "cProtein.hpp"


double scHemisphereRadii(cAminoResidue::eAminoType aminoType);

bool HaveCovalentBond(const cAtom &first, const cAtom &second);

std::map<cAtom *, std::set<cAtom *>> FindTopologicalNeighbours(cProtein *protein, size_t order);

std::map<cAtom *, std::set<cAtom *>> FindIsolatedNeighbours(cProtein *protein,
                                                            double cutoff,
                                                            size_t order);

template <typename T>
std::set<T> SetDifference(const std::set<T> &first, const std::set<T> &second) {
  std::set<T> result;
  std::set_difference(first.begin(), first.end(),
                      second.begin(), second.end(),
                      std::inserter(result, result.end()));
  return result;
}

std::map<cResidue *, std::vector<cResidue *>> NeighbourResidues(cProtein *protein,
                                                                double sidechain_cutoff,
                                                                double residue_dist_multiplier);

double contact_dist(cAminoResidue::eAminoType first, cAminoResidue::eAminoType second);

template <class Callback>
void ForEachAmino(cProtein *protein, Callback callback) {
  for (auto &residue : protein->residues()) {
    auto *res = dynamic_cast<cAminoResidue *>(&residue);
    if (!res)
      continue;

    callback(*res);
  }
}
