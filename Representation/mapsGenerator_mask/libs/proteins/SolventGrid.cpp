#include "SolventGrid.hpp"

#include <cstring>

#include "cGrid.hpp"
#include "cAlgorithmTemplates.hpp"
#include "cProtein.hpp"
#include "NeighbourSearch.hpp"


const double kHydrationShellResolution = 3.0;  // Angstrom
const double kPadding = 0.0;
const size_t kDummyWaterResidueSeqNumber = 999;
const char kDummyWaterResidueName[] = "HOH";
const size_t kDummyWaterOxygenSeqNumber = 999;


std::vector<cAtom *> GetOxygenGrid(cProtein *protein, double shell_width) {
  std::vector<cAtom *> atoms;
  for (cAtom &atom : protein->atoms()) {
    atoms.push_back(&atom);
  }

  double resolution = std::max(kHydrationShellResolution, shell_width);

  // resolution/times gives the spacing of the final Carteian grid
  int times = static_cast<int>(resolution / kHydrationShellResolution + 0.5);
  assert(times >= 1);

  //compute AABB around the protein
  cIAVector3 boundingBox = algorithms::computeAABB(atoms.data(), atoms.size());

  //expand the AABB by the kPadding distance
  if (kPadding > 0) {
    boundingBox.expand(kPadding);
  }
  // create a grid
  cGrid grid;

  // initialize the structure
  grid.initPointer(atoms.data(), atoms.size(), &boundingBox, resolution);
  grid.initLeafLayer(times); // inits leafs for non-occupied cells

  // first, basic version:

  // points that are not closer to protein atoms than VdW+solvent accessible radius
  // AND not further than resolution
  grid.detectBoundaryLeafs(0);

  cAtom water_pattern;
  water_pattern.atomSerial = kDummyWaterOxygenSeqNumber;
  std::strcpy(water_pattern.name, "O");
  water_pattern.altLoc = ' ';
  std::strcpy(water_pattern.resName, kDummyWaterResidueName);
  water_pattern.chainId = ' ';
  water_pattern.resSeqNumber = kDummyWaterResidueSeqNumber;
  water_pattern.insertion = ' ';
  std::strcpy(water_pattern.segId, "");
  std::strcpy(water_pattern.elementSymbol, "O");
  std::strcpy(water_pattern.chargeOnAtom, "");

  std::vector<cAtom *> water_atoms;
  water_atoms.reserve(grid.nVectors);

  for (int i = 0; i < grid.nVectors; i++) {
    cAtom *water_atom = new cAtom(water_pattern);
    water_atom->setPosition(*grid.listVectors[i]);
    water_atom->neighbours.clear();
    water_atoms.push_back(water_atom);
  }

  return water_atoms;
}

// Sets solvent weight to occupancy
void SetSolventOccupancy(cAtom *water_atom,
                         double sa_radius,
                         double sa_smoothed_margin) {
  // min dist from water to residue hemisphere
  double solvent_radius = sa_radius + sa_smoothed_margin;
  for (cAtom *neighbour : water_atom->neighbours) {
    double dist = (water_atom->getPosition() - neighbour->getPosition()).norm();
    auto *res = dynamic_cast<const cAminoResidue *>(neighbour->residue);
    if (dist - scHemisphereRadii(res->aminoType()) < solvent_radius) {
      solvent_radius = dist - scHemisphereRadii(res->aminoType());
    }
  }
  if (solvent_radius < sa_radius) {
    water_atom->occupancy = 0.0;
  } else if (solvent_radius > sa_radius + sa_smoothed_margin) {
    water_atom->occupancy = 1.0;
  } else {
    water_atom->occupancy = (solvent_radius - sa_radius) / sa_smoothed_margin;
  }
}

std::vector<std::unique_ptr<cWaterResidue>> GetHydrationShell(cProtein *protein,
                                                              double sa_radius,
                                                              double sa_smoothed_margin,
                                                              double shell_width) {
  auto oxygen_grid = GetOxygenGrid(protein, shell_width);

  // Set to zero occupancies for oxygen atoms that are too close to CA atoms
  std::vector<cAtom *> CA_atoms;
  ForEachAmino(protein, [&CA_atoms](cAminoResidue &residue) {
    auto *atom_CA = const_cast<cAtom *>(residue.getAtomByLabel(cAminoResidue::CA));
    assert(atom_CA);
    CA_atoms.push_back(atom_CA);
  });

  cGrid::initializeNeighbours(CA_atoms, oxygen_grid,
                              sa_radius + sa_smoothed_margin
                                        + scHemisphereRadii(cAminoResidue::ARG));

  std::vector<std::unique_ptr<cWaterResidue>> water_residues;

  for (cAtom *water_atom : oxygen_grid) {
    SetSolventOccupancy(water_atom, sa_radius, sa_smoothed_margin);
    if (water_atom->occupancy == 0) {
      delete water_atom;
      continue;
    }
    auto residue = cResidue::create(kDummyWaterResidueName, ' ', "",
                                    kDummyWaterResidueSeqNumber, NULL);
    cWaterResidue *water_residue = dynamic_cast<cWaterResidue *>(residue);
    assert(water_residue);

    water_residue->addAtom(water_atom);
    water_residues.push_back(std::unique_ptr<cWaterResidue>(water_residue));
  }
  return water_residues;
}
