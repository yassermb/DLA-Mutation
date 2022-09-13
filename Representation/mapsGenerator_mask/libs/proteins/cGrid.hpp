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


#pragma once

#include "cMemoryPool.hpp"
#include "cIAVector3.hpp"

#include <vector>
//struct sVoxel;
#include "cAtom.hpp"
#include "cProtein.hpp"

class cCell {
 public:
  virtual int getOccupied() = 0;
};


class  cCellLeafs: public cCell {
  friend class cGrid;
  // int occupied;
  int partitions;
  int partitions1D;
  cVector3 **vectors;
  virtual int getOccupied() {return -1;}

 public:
//    cCellLeafs() { contents = NULL; occupied=0;}
//    virtual ~cCellLeafs() {
//        if (contents) {
//            for (int i=0; i< partitions; i++)
//                delete contents[i];
//            delete [] contents;
//        }
//    }
};

class  cCellPDB: public cCellLeafs {
    friend class cGrid;
  int occupied;
  cAtom **contents;
  virtual int getOccupied() {return occupied;}
};

class cGrid {
 public:
  cGrid();
  ~cGrid();

  // Initialize neighbors for |atoms| among |can_be_neighbors|
  static void initializeNeighbours(std::vector<cAtom *> &can_be_neighbors,
                                   std::vector<cAtom *> &atoms,
                                   double cutoffDistance);

  // initialize the grid from an array of atoms
  void      initPointer(cAtom **_atoms, int _nAtoms, cIAVector3 *bb, double _cutoffDistance);

  int       contains(const cVector3 &pos);
  void      query(cAtom **inAtoms, int nIn, cAtom **outAtoms, int &nOut);
  void      queryOut(cAtom **inAtoms, int nIn, cAtom **outAtoms, int &nOut);
  void      constructNeighbourList(cAtom **array, int n);
  void		constructNeighbours(cVector3 position, std::set<cAtom*> &neighbours);
  void      initLeafLayer(int times);
  void      detectBoundaryLeafs(double cutoff);
  cVector3  **listVectors;
  int       nVectors;

 private:
  cMemoryPool pool;
  char *gridBit;
  cCell **gridPDB;

  cAtom **PDBatoms;
  int nAtoms;
  double cutoffDistance, inverseCutoff;
  cIAVector3 bound;
  int zyCells, zCells, xCells, yCells;
  int gridSize;
};
