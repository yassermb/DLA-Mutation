/*************************************************************************\

  Copyright 2004 The University of North Carolina at Chapel Hill.
  All Rights Reserved.

  Permission to use, copy, modify OR distribute this software and its
  documentation for educational, research and non-profit purposes, without
  fee, and without a written agreement is hereby granted, provided that the
  above copyright notice and the following three paragraphs appear in all
  copies.

  IN NO EVENT SHALL THE UNIVERSITY OF NORTH CAROLINA AT CHAPEL HILL BE
  LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR
  CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE
  USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF THE UNIVERSITY
  OF NORTH CAROLINA HAVE BEEN ADVISED OF THE POSSIBILITY OF SUCH
  DAMAGES.

  THE UNIVERSITY OF NORTH CAROLINA SPECIFICALLY DISCLAIM ANY
  WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.  THE SOFTWARE
  PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE UNIVERSITY OF
  NORTH CAROLINA HAS NO OBLIGATIONS TO PROVIDE MAINTENANCE, SUPPORT,
  UPDATES, ENHANCEMENTS, OR MODIFICATIONS.

  The author may be contacted via:

  US Mail:             Stephane Redon
                       Department of Computer Science
                       Sitterson Hall, CB #3175
                       University of N. Carolina
                       Chapel Hill, NC 27599-3175

  Phone:               (919) 962-1930

  EMail:               redon@cs.unc.edu

\**************************************************************************/


#pragma once

#include <cmath>

#include "cInterval.hpp"
#include "cVector3.hpp"

// implements a 3-dimensional interval vector

class cIAVector3 {

public:

	cInterval i[3]; // one interval per coordinate

	cIAVector3();
	cIAVector3(double val);
	cIAVector3(double x, double y, double z);
	cIAVector3(double xl, double xu, double yl, double yu, double zl, double zu);
	cIAVector3(cInterval val[3]);
	cIAVector3(double val[3][2]);
	cIAVector3(cInterval &nx, cInterval &ny, cInterval &nz);
	cIAVector3(cVector3 &u);
	cIAVector3(const cVector3 &u);

	cIAVector3		operator+(const cIAVector3 &u);
	cIAVector3&		operator+=(const cIAVector3 &u);
	cIAVector3&		operator+=(const cVector3 &v);
	cIAVector3		operator-(const cIAVector3 &u);
	cIAVector3&		operator-=(const cIAVector3 &u);
	cIAVector3&		operator=(const cVector3 &u);
	cInterval		operator|(const cIAVector3 &u);
	cIAVector3		operator^(const cIAVector3 &u);

	void			print();
	cVector3		center();
	double			volume();
	void			setZero();
	double			lengthsSum();
	void			makePointBox();
	cVector3		realComponent();
	void			contain(cIAVector3 &u);
	bool			overlaps(cIAVector3 &u);
	bool 			overlaps(cIAVector3 &u, double &cutoffDistance);
	bool			inBounds(cIAVector3 &u);
	double			distance2ToPoint(cVector3 &u);
	void			split(cIAVector3 &u1, cIAVector3 &u2);
	void			expand (double r);
	void			shrink (double r);

	enum 			{INCLUDE_ALL, INCLUDE_NONE, SUBDIVIDE};	
	unsigned int 	overlapsAdvanced(cIAVector3 &u, const double &cutoff, const double &cutoff2); 

	void			bound(const cIAVector3 &u);
	
};

typedef cIAVector3 *pcIAVector3;


