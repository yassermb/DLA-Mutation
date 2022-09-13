/*************************************************************************\

 Copyright 2010 Sergei Grudinin, INRIA.
 All Rights Reserved.

 The author may be contacted via:

 Mail:          Sergei Grudinin
 INRIA Rhone-Alpes Research Unit
 Zirst - 655 avenue de l'Europe - Montbonnot
 38334 Saint Ismier Cedex - France

 Phone:       +33 4 76 61 53 24

 EMail:       sergei.grudinin@inria.fr

 \**************************************************************************/


#include <string.h>
#include <assert.h>

#include "cGrid.hpp"
#include "cMemoryPool.hpp"
#include "cIAVector3.hpp"
#include "cAlgorithmTemplates.hpp"
#include "cAtom.hpp"
#include "cProtein.hpp"

cGrid::~cGrid() {

  if (gridBit)
    delete [] gridBit;
  if (gridPDB)
    delete [] gridPDB;
  if (listVectors)
    delete [] listVectors;
}

cGrid::cGrid() : pool(1024*1024) {

  gridBit = NULL;
  gridPDB = NULL;

  listVectors = NULL;
  nVectors = 0;
}

void cGrid::initializeNeighbours(std::vector<cAtom *> &can_be_neighbours,
                                 std::vector<cAtom *> &atoms,
                                 double cutoffDistance) {
  if (!atoms.size() || !can_be_neighbours.size())
    return;

  cGrid grid;
  grid.initPointer(can_be_neighbours.data(),
                   can_be_neighbours.size(),
                   NULL,
                   cutoffDistance);
  grid.constructNeighbourList(atoms.data(), atoms.size());
}

int cGrid::contains(const cVector3 &pos) {

  unsigned int x,y,z,ind, bitInd;
  char mask;
  x = (pos[0]-bound.i[0].i[0])*inverseCutoff+1;
  y = (pos[1]-bound.i[1].i[0])*inverseCutoff+1;
  z = (pos[2]-bound.i[2].i[0])*inverseCutoff+1;

  ind = (x*zyCells) + (y*zCells) + z;
  bitInd = ind >> 3;
  mask = 1<<(ind&7);

  return gridBit[bitInd]&mask;
}

void cGrid::query(cAtom **inAtoms, int nIn, cAtom **outAtoms, int &nOut) {

  for (int i=0; i<nIn; i++) {
    if (this->contains(inAtoms[i]->getPosition())) {
      outAtoms[nOut++] = inAtoms[i];
    }
  }
}

void cGrid::queryOut(cAtom **inAtoms, int nIn, cAtom **outAtoms, int &nOut) {

  for (int i=0; i<nIn; i++) {
    if(!this->contains(inAtoms[i]->getPosition())) {
      outAtoms[nOut++] = inAtoms[i];
    }
  }
}


#define COMPARE_FUNC_PDB(k) \
currentCell = static_cast<cCellPDB*>(gridPDB[ind[k]]); \
if (currentCell) { \
inside = currentCell->occupied; \
for (int j=0; j<=inside; j++) { \
const cVector3  &tmpV = currentCell->contents[j]->getPosition(); \
if ((curPosition - tmpV).norm2() < cutoff2 ) { \
neighbours.insert(currentCell->contents[j]); \
} \
} \
}

void cGrid::constructNeighbours(cVector3 curPosition, std::set<cAtom*> &neighbours) {
	// +1 for cutoffl
	// +1 for padding
	double cutoff2 = cutoffDistance*cutoffDistance;
	cCellPDB* currentCell;
	int x = (curPosition[0] - bound.i[0].i[0])*inverseCutoff + 2;
	int y = (curPosition[1] - bound.i[1].i[0])*inverseCutoff + 2;
	int z = (curPosition[2] - bound.i[2].i[0])*inverseCutoff + 2;
	int inside;

	// outside the bounding box
	if (x<1 || x >= (xCells - 1) || y<1 || y >= (yCells - 1) || z<1 || z >= (zCells - 1))
		return;

	int index = (x*zyCells) + (y*zCells) + z;

	// look around and inside

	int ind[27];
	ind[0] = index - 1;
	ind[1] = index + 1;
	ind[2] = index + zCells;
	ind[3] = index - zCells;
	ind[4] = ind[2] + 1;
	ind[5] = ind[2] - 1;
	ind[6] = ind[3] + 1;
	ind[7] = ind[3] - 1;
	//
	ind[8] = index + zyCells;
	ind[9] = ind[8] + 1;
	ind[10] = ind[8] - 1;
	ind[11] = ind[8] + zCells;
	ind[12] = ind[8] - zCells;
	ind[13] = ind[11] + 1;
	ind[14] = ind[11] - 1;
	ind[15] = ind[12] + 1;
	ind[16] = ind[12] - 1;
	//
	ind[17] = index - zyCells;
	ind[18] = ind[17] + 1;
	ind[19] = ind[17] - 1;
	ind[20] = ind[17] + zCells;
	ind[21] = ind[17] - zCells;
	ind[22] = ind[20] + 1;
	ind[23] = ind[20] - 1;
	ind[24] = ind[21] + 1;
	ind[25] = ind[21] - 1;
	//
	ind[26] = index;


	// unroll a bit...
	COMPARE_FUNC_PDB(0);
	COMPARE_FUNC_PDB(1);
	COMPARE_FUNC_PDB(2);
	COMPARE_FUNC_PDB(3);
	COMPARE_FUNC_PDB(4);
	COMPARE_FUNC_PDB(5);
	COMPARE_FUNC_PDB(6);
	COMPARE_FUNC_PDB(7);
	COMPARE_FUNC_PDB(8);
	COMPARE_FUNC_PDB(9);
	COMPARE_FUNC_PDB(10);
	COMPARE_FUNC_PDB(11);
	COMPARE_FUNC_PDB(12);
	COMPARE_FUNC_PDB(13);
	COMPARE_FUNC_PDB(14);
	COMPARE_FUNC_PDB(15);
	COMPARE_FUNC_PDB(16);
	COMPARE_FUNC_PDB(17);
	COMPARE_FUNC_PDB(18);
	COMPARE_FUNC_PDB(19);
	COMPARE_FUNC_PDB(20);
	COMPARE_FUNC_PDB(21);
	COMPARE_FUNC_PDB(22);
	COMPARE_FUNC_PDB(23);
	COMPARE_FUNC_PDB(24);
	COMPARE_FUNC_PDB(25);
	COMPARE_FUNC_PDB(26);

}

void cGrid::constructNeighbourList(cAtom **array, int n) {

  cVector3 curPosition;
  int x,y,z, inside, index;
  double cutoff2 = cutoffDistance*cutoffDistance;
  cCellPDB* currentCell;

  for (int i = 0; i < n; i++) {

    curPosition= array[i]->getPosition();
	constructNeighbours(curPosition, array[i]->neighbours);
	array[i]->neighbours.erase(array[i]);

  }
}

// initialize the grid from an array of atoms
void cGrid::initPointer(cAtom **_atoms, int _nAtoms, cIAVector3 *bb, double _cutoffDistance) {
	assert(_nAtoms > 0);

	PDBatoms = _atoms;
	nAtoms = _nAtoms;
	cutoffDistance = _cutoffDistance;

	cVector3 curPosition;
	inverseCutoff = 1.0 / cutoffDistance;

	// compute the Bounding Box if empty
	if (!bb) {
		bound = algorithms::computeAABB(PDBatoms, nAtoms);
	}
	else {
		bound = *bb;
	}

	// +3 instead of +1 (each side adds cutoff, don't forget it when will put)
	// (int) is much faster than floor ~ 4-50 times
	// + additional 2 to prevent if-stuff
	xCells = (int)((bound.i[0].diameter())*inverseCutoff) + 3 + 2;
	yCells = (int)((bound.i[1].diameter())*inverseCutoff) + 3 + 2;
	zCells = (int)((bound.i[2].diameter())*inverseCutoff) + 3 + 2;

	zyCells = yCells*zCells;
	gridSize = xCells*zyCells;

	gridPDB = new cCell*[gridSize + 1]; //+1 to have NULL at the end
	memset(gridPDB, 0, sizeof(cCellPDB*)*(gridSize + 1));

	// how much can be inside, probably another strategy is better
	int inside = (3 * cutoffDistance*cutoffDistance*cutoffDistance);

	unsigned int x, y, z, ind;

	for (int i = 0; i<nAtoms; i++) {

		curPosition = PDBatoms[i]->getPosition();

		//+ how many?
		// +1 for cutoff
		// +1 for padding
		x = (curPosition[0] - bound.i[0].i[0])*inverseCutoff + 2;
		y = (curPosition[1] - bound.i[1].i[0])*inverseCutoff + 2;
		z = (curPosition[2] - bound.i[2].i[0])*inverseCutoff + 2;

		ind = (x*zyCells) + (y*zCells) + z;

		//cCellPDB* currentCell = static_cast<cCellPDB*>(gridPDB[ind]);
		//cCell* currentCell = gridPDB[ind];


		if (gridPDB[ind] != NULL) {

			cCellPDB* currentCell = static_cast<cCellPDB*>(gridPDB[ind]);
			currentCell->occupied++;
			currentCell->contents[currentCell->occupied] = PDBatoms[i];

			assert(currentCell->occupied < inside);

		}
		else {

			gridPDB[ind] = new(pool) cCellPDB;
			cCellPDB* currentCell = static_cast<cCellPDB*>(gridPDB[ind]);
			currentCell->contents = new(pool) cAtom*[inside];
			currentCell->contents[0] = PDBatoms[i];
			currentCell->occupied = 0;
		}

	}

}



void cGrid::initLeafLayer(int times) {

    int times2 = times*times;
    int times3 = times2*times;

    listVectors = new cVector3*[(xCells-2)*(yCells-2)*(zCells-2)*times3];

    for (int x = 1; x < xCells-1; x++)
        for (int y = 1; y < yCells-1; y++)
            for (int z = 1; z < zCells-1; z++) {
               int index = (x*zyCells) + (y*zCells) + z;
                cCell *&cell = gridPDB[index];
                if (!cell) {
                    cCellLeafs * newCell = new(pool) cCellLeafs();
                    cell = newCell;
                }
                cCellLeafs * newCell = static_cast<cCellLeafs *>(cell);
                //    x = (curPosition[0]-bound.i[0].i[0])*inverseCutoff+2;
                double cellX = (x-2)*cutoffDistance + bound.i[0].i[0];
                double cellY = (y-2)*cutoffDistance + bound.i[1].i[0];
                double cellZ = (z-2)*cutoffDistance + bound.i[2].i[0];

                newCell->partitions = times3;
                newCell->partitions1D = times;
                newCell->vectors = new(pool) cVector3*[newCell->partitions];
                for ( int i1=0; i1< times; i1++) {
                    double newX = cellX + cutoffDistance/times*(i1+0.5);
                    for ( int i2=0; i2< times; i2++) {
                        double newY = cellY + cutoffDistance/times*(i2+0.5);
                        for ( int i3=0; i3< times; i3++) {
                            double newZ = cellZ + cutoffDistance/times*(i3+0.5);
                            int ind = i1*times2+ i2*times + i3;
                            newCell->vectors[ind] = new(pool) cVector3(newX,newY,newZ);
                        }
                    }
                }

            }


}

#define COMPARE_FUNC_LEAF(k) \
cell =  gridPDB[ind[k]]; \
if (cell && (cell->getOccupied() != -1)) { \
cellPDB = static_cast<cCellPDB*>(cell); \
inside = cellPDB->occupied; \
for (int j=0; j<=inside; j++) { \
float R = cellPDB->contents[j]->getElementRadius(); \
cutoffPDB2 = (cutoff+R) *(cutoff+R);\
const cVector3  &tmpV = cellPDB->contents[j]->getPosition(); \
if (( distance2 = (*curPosition - tmpV).norm2()) < cutoffPDB2) { \
goto loop; \
} else if (distance2 < cutoff2) {\
++insideCutoffDistance; \
} \
} \
}

void cGrid::detectBoundaryLeafs(double cutoff) {

  cCell *cell;
  cCellPDB *cellPDB;
  int inside;
  double cutoffPDB2, distance2;
  double cutoff2 = cutoffDistance*cutoffDistance;
  // int k;

  for (int x = 1; x < xCells-1; x++) {
    for (int y = 1; y < yCells-1; y++) {
      for (int z = 1; z < zCells-1; z++) {
        int index = (x*zyCells) + (y*zCells) + z;
        cCell *currentCell = gridPDB[index];
        // if it's not empty:
        // if (currentCell->getOccupied() != -1)
          // continue;

        // look around
        int ind[27];
        ind[0]=index-1;
        ind[1]=index+1;
        ind[2]=index+zCells;
        ind[3]=index-zCells;
        ind[4]=ind[2]+1;
        ind[5]=ind[2]-1;
        ind[6]=ind[3]+1;
        ind[7]=ind[3]-1;
        //
        ind[8]=index+zyCells;
        ind[9]=ind[8]+1;
        ind[10]=ind[8]-1;
        ind[11]=ind[8]+zCells;
        ind[12]=ind[8]-zCells;
        ind[13]=ind[11]+1;
        ind[14]=ind[11]-1;
        ind[15]=ind[12]+1;
        ind[16]=ind[12]-1;
        //
        ind[17]=index-zyCells;
        ind[18]=ind[17]+1;
        ind[19]=ind[17]-1;
        ind[20]=ind[17]+zCells;
        ind[21]=ind[17]-zCells;
        ind[22]=ind[20]+1;
        ind[23]=ind[20]-1;
        ind[24]=ind[21]+1;
        ind[25]=ind[21]-1;
        ind[26]=index;


        int times = static_cast<cCellLeafs*>(currentCell)->partitions1D;
        int times2 = times*times;
        for ( int i1=0; i1< times; i1++) {
          for ( int i2=0; i2< times; i2++) {
            for ( int i3=0; i3< times; i3++) {
              int indCell = i1*times2+ i2*times + i3;
              cVector3 *curPosition
                = static_cast<cCellLeafs*>(currentCell)->vectors[indCell];

              int insideCutoffDistance = 0;
              COMPARE_FUNC_LEAF(0);
              COMPARE_FUNC_LEAF(1);
              COMPARE_FUNC_LEAF(2);
              COMPARE_FUNC_LEAF(3);
              COMPARE_FUNC_LEAF(4);
              COMPARE_FUNC_LEAF(5);
              COMPARE_FUNC_LEAF(6);
              COMPARE_FUNC_LEAF(7);
              COMPARE_FUNC_LEAF(8);
              COMPARE_FUNC_LEAF(9);
              COMPARE_FUNC_LEAF(10);
              COMPARE_FUNC_LEAF(11);
              COMPARE_FUNC_LEAF(12);
              COMPARE_FUNC_LEAF(13);
              COMPARE_FUNC_LEAF(14);
              COMPARE_FUNC_LEAF(15);
              COMPARE_FUNC_LEAF(16);
              COMPARE_FUNC_LEAF(17);
              COMPARE_FUNC_LEAF(18);
              COMPARE_FUNC_LEAF(19);
              COMPARE_FUNC_LEAF(20);
              COMPARE_FUNC_LEAF(21);
              COMPARE_FUNC_LEAF(22);
              COMPARE_FUNC_LEAF(23);
              COMPARE_FUNC_LEAF(24);
              COMPARE_FUNC_LEAF(25);
//DEBUG only
              // k=25;
              // cell =  gridPDB[ind[k]];
              // if (cell && (cell->getOccupied() >= 0)) {
                 // cellPDB = static_cast<cCellPDB*>(cell);
                 // inside = cellPDB->occupied;
                 // for (int j=0; j<=inside; j++) {
                     // float &R = cellPDB->contents[j]->elementRadius;
                     // cutoffPDB2 = (cutoff+R) *(cutoff+R);
                     // const cVector3  &tmpV = cellPDB->contents[j]->getPosition();
                     // if (( ((curPosition->v[0] - tmpV[0])*(curPosition->v[0] - tmpV[0]) +(curPosition->v[1] - tmpV[1])*(curPosition->v[1] - tmpV[1]) +(curPosition->v[2] - tmpV[2])*(curPosition->v[2] - tmpV[2]))) < cutoffPDB2) {
                         // goto loop;
                     // }
                 // }
              // }
              COMPARE_FUNC_LEAF(26);

              if (insideCutoffDistance)
                listVectors[nVectors++]=curPosition;
              loop:
              ;
            }
          }
        }
      }
    }
  }
}
