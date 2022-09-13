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


#include <math.h>
#include <cstring>

#include "cIAVector3.hpp"
#include "cInterval.hpp"
#include "cVector3.hpp"


#ifndef signbit
#define signbit(A) ( (A)<0 ? true : false )
#endif

cIAVector3::cIAVector3() {}
cIAVector3::cIAVector3(double val) { i[0]=i[1]=i[2]=val; }
cIAVector3::cIAVector3(double x, double y, double z) { i[0].i[0]=x;i[0].i[0]=y;i[0].i[0]=z; }
cIAVector3::cIAVector3(double xl, double xu, double yl, double yu, double zl, double zu) {
	i[0].i[0]=xl;i[0].i[1]=xu;
	i[1].i[0]=yl;i[1].i[1]=yu;
	i[2].i[0]=zl;i[2].i[1]=zu;
}
cIAVector3::cIAVector3(double val[3][2]) { memcpy(i,val,6*sizeof(double)); }
cIAVector3::cIAVector3(cInterval val[3]) { memcpy(i,val,6*sizeof(double)); }
cIAVector3::cIAVector3(cInterval &nx, cInterval &ny, cInterval &nz) { i[0]=nx;i[1]=ny;i[2]=nz; }
cIAVector3::cIAVector3(cVector3 &u) { i[0]=u[0];i[1]=u[1];i[2]=u[2]; }
cIAVector3::cIAVector3(const cVector3 &u) { i[0]=u[0];i[1]=u[1];i[2]=u[2]; }

void cIAVector3::setZero() { i[0]=i[1]=i[2]=0.0; }

// operators

cIAVector3 cIAVector3::operator+(const cIAVector3 &u) { cInterval res[3];res[0]=i[0]+u.i[0];res[1]=i[1]+u.i[1];res[2]=i[2]+u.i[2];return cIAVector3(res); }
cIAVector3& cIAVector3::operator+=(const cIAVector3 &u) { i[0]+=u.i[0];i[1]+=u.i[1];i[2]+=u.i[2];return *this; }
cIAVector3 cIAVector3::operator-(const cIAVector3 &u) { cInterval res[3];res[0]=i[0]-u.i[0];res[1]=i[1]-u.i[1];res[2]=i[2]-u.i[2];return cIAVector3(res); }
cIAVector3& cIAVector3::operator-=(const cIAVector3 &u) { i[0]-=u.i[0];i[1]-=u.i[1];i[2]-=u.i[2];return *this; }
cIAVector3& cIAVector3::operator=(const cVector3 &u) { i[0]=u[0];i[1]=u[1];i[2]=u[2];return *this; }
cInterval cIAVector3::operator|(const cIAVector3 &u) { return i[0]*u.i[0]+i[1]*u.i[1]+i[2]*u.i[2]; }
cIAVector3 cIAVector3::operator^(const cIAVector3 &u) { cInterval res[3];res[0]=i[1]*u.i[2]-i[2]*u.i[1];res[1]=i[2]*u.i[0]-i[0]*u.i[2];res[2]=i[0]*u.i[1]-i[1]*u.i[0];return cIAVector3(res); }
cIAVector3& cIAVector3::operator+=(const cVector3 &v) { i[0]+=v[0];i[1]+=v[1];i[2]+=v[2];return *this; }
double cIAVector3::lengthsSum() { return i[0].diameter()+i[1].diameter()+i[2].diameter(); }
double cIAVector3::volume() { return i[0].diameter()*i[1].diameter()*i[2].diameter(); }
void cIAVector3::makePointBox() { i[0].i[1]=i[0].i[0];i[1].i[1]=i[1].i[0];i[2].i[1]=i[2].i[0]; }
void cIAVector3::print() {
	std::cout << " [" << i[0].i[0] << "," << i[0].i[1] << "]" << std::endl;
	std::cout << " [" << i[1].i[0] << "," << i[1].i[1] << "]" << std::endl;
	std::cout << " [" << i[2].i[0] << "," << i[2].i[1] << "]" << std::endl;
	std::cout << std::endl;
}
cVector3 cIAVector3::center() { return cVector3(i[0].center(),i[1].center(),i[2].center()); }
void cIAVector3::contain(cIAVector3 &u) { // enlarge the AABB represented by the IA std::vector to contain the one represented by u
	
	if (u.i[0].i[0]<i[0].i[0]) i[0].i[0]=u.i[0].i[0];
	if (u.i[1].i[0]<i[1].i[0]) i[1].i[0]=u.i[1].i[0];
	if (u.i[2].i[0]<i[2].i[0]) i[2].i[0]=u.i[2].i[0];
	if (u.i[0].i[1]>i[0].i[1]) i[0].i[1]=u.i[0].i[1];
	if (u.i[1].i[1]>i[1].i[1]) i[1].i[1]=u.i[1].i[1];
	if (u.i[2].i[1]>i[2].i[1]) i[2].i[1]=u.i[2].i[1];
	
}
bool cIAVector3::overlaps(cIAVector3 &u) { // returns true if and only if the AABB u overlaps the AABB
	
	if (u.i[0].i[1]<i[0].i[0]) return false;
	if (u.i[0].i[0]>i[0].i[1]) return false;
	if (u.i[1].i[1]<i[1].i[0]) return false;
	if (u.i[1].i[0]>i[1].i[1]) return false;
	if (u.i[2].i[1]<i[2].i[0]) return false;
	if (u.i[2].i[0]>i[2].i[1]) return false;
	
	return true;
	
}

bool cIAVector3::overlaps(cIAVector3 &u, double &cutoffDistance) { // returns true if and only if the AABB u overlaps the AABB
	
	if (u.i[0].i[1]+cutoffDistance<i[0].i[0]) return false;
	if (u.i[0].i[0]>i[0].i[1]+cutoffDistance) return false;
	if (u.i[1].i[1]+cutoffDistance<i[1].i[0]) return false;
	if (u.i[1].i[0]>i[1].i[1]+cutoffDistance) return false;
	if (u.i[2].i[1]+cutoffDistance<i[2].i[0]) return false;
	if (u.i[2].i[0]>i[2].i[1]+cutoffDistance) return false;
	
	return true;
	
}

cVector3 cIAVector3::realComponent() { return cVector3(i[0].i[0],i[1].i[0],i[2].i[0]); }

bool cIAVector3::inBounds(cIAVector3 &u) {

	if (i[0].i[0]<u.i[0].i[0]) return false;
	if (i[0].i[1]>u.i[0].i[1]) return false;
	if (i[1].i[0]<u.i[1].i[0]) return false;
	if (i[1].i[1]>u.i[1].i[1]) return false;
	if (i[2].i[0]<u.i[2].i[0]) return false;
	if (i[2].i[1]>u.i[2].i[1]) return false;

	return true;

}

double cIAVector3::distance2ToPoint(cVector3 &u) {

	double x,y,z;

	if (u[0]<=i[0].i[0]) x=i[0].i[0]; 
	else if (u[0]>=i[0].i[1]) x=i[0].i[1]; 
	else x=u[0];

	if (u[1]<=i[1].i[0]) y=i[1].i[0]; 
	else if (u[1]>=i[1].i[1]) y=i[1].i[1]; 
	else y=u[1];

	if (u[2]<=i[2].i[0]) z=i[2].i[0]; 
	else if (u[2]>=i[2].i[1]) z=i[2].i[1]; 
	else z=u[2];

	return (x-u[0])*(x-u[0])+(y-u[1])*(y-u[1])+(z-u[2])*(z-u[2]);

}

void cIAVector3::split(cIAVector3 &u1, cIAVector3 &u2) {

	// split along the longest axis

	int axis=0;
	double maxLength=i[0].diameter();
	double m;

	if (i[1].diameter()>maxLength) { axis=1;maxLength=i[1].diameter(); }
	if (i[2].diameter()>maxLength) { axis=2;maxLength=i[2].diameter(); }

	switch (axis) {

		case 0:

			m=i[0].center();
			u1.i[0].i[0]=i[0].i[0];u1.i[0].i[1]=m;
			u2.i[0].i[1]=i[0].i[1];u2.i[0].i[0]=m;
			u1.i[1]=i[1];u2.i[1]=i[1];
			u1.i[2]=i[2];u2.i[2]=i[2];
			break;

		case 1:

			m=i[1].center();
			u1.i[1].i[0]=i[1].i[0];u1.i[1].i[1]=m;
			u2.i[1].i[1]=i[1].i[1];u2.i[1].i[0]=m;
			u1.i[0]=i[0];u2.i[0]=i[0];
			u1.i[2]=i[2];u2.i[2]=i[2];
			break;

		case 2:

			m=i[2].center();
			u1.i[2].i[0]=i[2].i[0];u1.i[2].i[1]=m;
			u2.i[2].i[1]=i[2].i[1];u2.i[2].i[0]=m;
			u1.i[0]=i[0];u2.i[0]=i[0];
			u1.i[1]=i[1];u2.i[1]=i[1];
			break;

	}

}

void cIAVector3::expand (double r) {

	i[0].expand(r);
	i[1].expand(r);
	i[2].expand(r);
}

void cIAVector3::shrink (double r) {
	
	i[0].shrink(r);
	i[1].shrink(r);
	i[2].shrink(r);
}


unsigned int cIAVector3::overlapsAdvanced(cIAVector3 &u, const double &cutoff, const double &cutoff2) { // returns true if and only if the AABB u overlaps the AABB

	bool signA, signB;
	double max2=0.0;
	double min2=0.0;
	double a, b;
	
	if ((a = u.i[0].i[1] - i[0].i[0]) < -cutoff)
		return INCLUDE_NONE; 
	if ((b = u.i[0].i[0] - i[0].i[1]) > cutoff)
		return INCLUDE_NONE; 
	
	if ((signA = signbit(a))^(signB = signbit(b)))  { // intersection, sign(a0)=true, sign(b0)=false
		max2 += a>fabs(b) ? a*a: b*b;
		//min2 += 0;
	}
	else if (!signA) {
			max2 += a*a; min2 += b*b;
	} 
	else {
			max2 += b*b; min2 += a*a;
	}
	
	if ((a = u.i[1].i[1] - i[1].i[0]) < -cutoff)
		return INCLUDE_NONE; 
	if ((b = u.i[1].i[0] - i[1].i[1]) > cutoff)
		return INCLUDE_NONE; 

	
	if ((signA = signbit(a))^(signB = signbit(b)))  { // intersection, sign(a0)=true, sign(b0)=false
		max2 += a>fabs(b) ? a*a: b*b;
		//min2 += 0;
	}
	else if (!signA) {
			max2 += a*a; min2 += b*b;
	} 
	else {
			max2 += b*b; min2 += a*a;
	}
		
	if ((a = u.i[2].i[1] - i[2].i[0]) < -cutoff)
		return INCLUDE_NONE; 
	if ((b = u.i[2].i[0] - i[2].i[1]) > cutoff)
		return INCLUDE_NONE; 

	
	if ((signA = signbit(a))^(signB = signbit(b)))  { // intersection, sign(a0)=true, sign(b0)=false
		max2 += a>fabs(b) ? a*a: b*b;
		//min2 += 0;
	}
	else if (!signA) {
			max2 += a*a; min2 += b*b;
	} 
	else {
			max2 += b*b; min2 += a*a;
	}
		
	if (max2 < cutoff2) // include all pairs
		return  INCLUDE_ALL;
	else if (min2 > cutoff2) // none is included
		return INCLUDE_NONE; 
	else	// subdivide
		return SUBDIVIDE;
	
}

void			cIAVector3::bound(const cIAVector3 &u) {

	i[0].bound(u.i[0]);
	i[1].bound(u.i[1]);
	i[2].bound(u.i[2]);

}
