/*************************************************************************\

 Copyright 2010 Sergei Grudinin, INRIA.
 All Rights Reserved.

 The author may be contacted via:

 Mail:					Sergei Grudinin
 INRIA Rhone-Alpes Research Unit
 Zirst - 655 avenue de l'Europe - Montbonnot
 38334 Saint Ismier Cedex - France

 Phone:				+33 4 76 61 53 24

 EMail:				sergei.grudinin@inria.fr

 \**************************************************************************/

#pragma once

#include <assert.h>
#include <algorithm>

#include "cIAVector3.hpp"


namespace algorithms {

	template <class dataType, class sizeType>
	void reset3D(dataType*** data, sizeType size) {

		for (unsigned int i = 0; i < size[0]; i++) {
			for (unsigned int j = 0; j < size[1]; j++) {
			 	memset(data[i][j], 0, size[2] * sizeof(dataType));
			}
		}

	}

	template <class dataType , class sizeType>
	void addArray3D(dataType*** dataOut, dataType*** data, sizeType size) {

		for (unsigned int i = 0; i < size[0]; i++) {
			for (unsigned int j = 0; j < size[1]; j++)
				for (unsigned int k = 0; k < size[2]; k++)
					dataOut[i][j][k] += data[i][j][k];

		}

	}

	template <class atomTemplate>
	bool compFunction_0(const atomTemplate &i, const atomTemplate &j) {
		return i.getPosition().v[0] < j.getPosition().v[0];
	}

	template <class atomTemplate>
	bool compFunction_1(const atomTemplate &i, const atomTemplate &j) {
		return i.getPosition().v[1] < j.getPosition().v[1];
	}

	template <class atomTemplate>
	bool compFunction_2(const atomTemplate &i, const atomTemplate &j) {
		return i.getPosition().v[2] < j.getPosition().v[2];
	}

	template <class atomTemplate>
	bool compFunction0(const atomTemplate *i, const atomTemplate *j) {
		return i->getPosition().v[0] < j->getPosition().v[0];
	}

	template <class atomTemplate>
	bool compFunction1(const atomTemplate *i, const atomTemplate *j) {
		return i->getPosition().v[1] < j->getPosition().v[1];
	}

	template <class atomTemplate>
	bool compFunction2(const atomTemplate *i, const atomTemplate *j) {
		return i->getPosition().v[2] < j->getPosition().v[2];
	}

	template <class containerTemplate, class atomTemplate>
	containerTemplate* createNonbondedAssemblyTree(atomTemplate **nodes, int n) {

		assert(n > 0);

		if (n == 1)
			return nodes[0]->getContainer();

		cVector3 center(0.0, 0.0, 0.0);

		if (n == 2) {
			center = (nodes[0]->getPosition() + nodes[1]->getPosition()) * 0.5;
		} else {

			double maxX, minX, maxY, minY, maxZ, minZ;
			maxX = minX = nodes[0]->getPosition().v[0];
			maxY = minY = nodes[0]->getPosition().v[1];
			maxZ = minZ = nodes[0]->getPosition().v[2];

			for (int i = 0; i < n; i++) {
				center += nodes[i]->getPosition();

				if (maxX < (nodes[i]->getPosition()).v[0])
					maxX = nodes[i]->getPosition().v[0];

				if (maxY < nodes[i]->getPosition().v[1])
					maxY = nodes[i]->getPosition().v[1];

				if (maxZ < nodes[i]->getPosition().v[2])
					maxZ = nodes[i]->getPosition().v[2];

				if (minX > nodes[i]->getPosition().v[0])
					minX = nodes[i]->getPosition().v[0];

				if (minY > nodes[i]->getPosition().v[1])
					minY = nodes[i]->getPosition().v[1];

				if (minZ > nodes[i]->getPosition().v[2])
					minZ = nodes[i]->getPosition().v[2];
			}

			int axis;

			if ((maxY - minY) > (maxX - minX)) {
				axis = 1;
				if ((maxZ - minZ) > (maxY - minY))
					axis = 2;
			}
			else {
				axis = 0;
				if ((maxZ - minZ) > (maxX - minX))
					axis = 2;
			}

			center /= static_cast<double>(n);

			switch (axis) {
				case 0:
					std::nth_element(nodes, nodes + n/2, nodes + n, compFunction0<atomTemplate>);
					break;
				case 1:
					std::nth_element(nodes, nodes + n/2, nodes + n, compFunction1<atomTemplate>);
					break;
				case 2:
					std::nth_element(nodes, nodes + n/2, nodes + n, compFunction2<atomTemplate>);
					break;
			}

		}

		containerTemplate* bodyA;
		containerTemplate* bodyB;

		atomTemplate **second = nodes + n/2;
		bodyA = createNonbondedAssemblyTree<containerTemplate, atomTemplate>(nodes, n/2);
		bodyB = createNonbondedAssemblyTree<containerTemplate, atomTemplate>(second, n - n/2);

		// gather together and return;

		containerTemplate* bodyC = new containerTemplate(bodyA, bodyB);

		// and assign the center of mass
		//bodyC->position = center;

		return bodyC;

	}

	template <class atomTempate>
	cIAVector3 computeAABB(atomTempate **atoms, int n) {
		assert(n > 0);

		cIAVector3 out(atoms[0]->getPosition());

		for (int i = 1; i < n; i++) {
			out.bound(cIAVector3(atoms[i]->getPosition()));
		}

		return out;
	}

	template <class atomTempate>
	cIAVector3 computeAABB(atomTempate *atoms, int n) {
		assert(n > 0);

		cIAVector3 out(atoms[0].getPosition());

		for (int i = 1; i < n; i++) {
			out.bound(cIAVector3(atoms[i].getPosition()));
		}

		return out;
	}

	template <class atomIt>
	cIAVector3 computeAABB(atomIt atoms, int n) {
		assert(n > 0);

		cIAVector3 out(atoms->getPosition());

		for (int i = 1; i < n; i++) {
			out.bound(cIAVector3((atoms+i)->getPosition()));
		}

		return out;
	}

	template <class atomIt>
	cIAVector3 computeAABB(atomIt &atom_begin, atomIt &atom_end) {
		assert(atom_begin != atom_end);

		cIAVector3 out(atom_begin->getPosition());

		for (atomIt it = atom_begin; it != atom_end; ++it) {
			out.bound(cIAVector3(it->getPosition()));
		}

		return out;
	}

	template <class atomTempate>
	void moveArrayBy(atomTempate **atoms, int n, cVector3 shift) {

		for (int i = 0; i < n; i++) {
			cVector3 coords = atoms[i]->getPosition();
			atoms[i]->setPosition(coords + shift);
		}

	}

	template <class atomTempate>
	void moveArrayBy(atomTempate *atoms, int n, cVector3 shift) {

		for (int i = 0; i < n; i++) {
			cVector3 coords = atoms[i].getPosition();
			atoms[i].setPosition(coords + shift);
		}

	}
}
