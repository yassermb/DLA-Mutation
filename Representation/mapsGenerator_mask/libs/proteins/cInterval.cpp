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

#include <iostream>

#include "cInterval.hpp"

cInterval::cInterval() {; }
cInterval::cInterval(double v) { i[0]=i[1]=v; }
cInterval::cInterval(double ll, double rr) { i[0]=ll;i[1]=rr; }
cInterval cInterval::operator+(const cInterval &in) { return cInterval(i[0]+in.i[0],i[1]+in.i[1]); }
cInterval& cInterval::operator+=(const double d) { i[0]+=d; i[1]+=d; return *this; }
cInterval cInterval::operator-(const cInterval &in) { return cInterval(i[0]-in.i[1],i[1]-in.i[0]); }
cInterval& cInterval::operator+=(const cInterval &in) { i[0]+=in.i[0];i[1]+=in.i[1];return *this; }
cInterval& cInterval::operator-=(const cInterval &in) { i[0]-=in.i[1];i[1]-=in.i[0];return *this; }
cInterval cInterval::operator*(const cInterval &in) { 
	
	if (in.i[0]>=0) {
		if (i[0]>=0) return cInterval(i[0]*in.i[0],i[1]*in.i[1]);
		if (i[1]<=0) return cInterval(i[0]*in.i[1],i[1]*in.i[0]);
		return cInterval(i[0]*in.i[1],i[1]*in.i[1]);
	}
	if (in.i[1]<=0) {
		if (i[0]>=0) return cInterval(i[1]*in.i[0],i[0]*in.i[1]);
		if (i[1]<=0) return cInterval(i[1]*in.i[1],i[0]*in.i[0]);
		return cInterval(i[1]*in.i[0],i[0]*in.i[0]);
	}
	if (i[0]>=0) return cInterval(i[1]*in.i[0],i[1]*in.i[1]);
	if (i[1]<=0) return cInterval(i[0]*in.i[1],i[0]*in.i[0]);
	double v00=i[0]*in.i[0];
	double v11=i[1]*in.i[1];
	if (v00<=v11) {
		double v01=i[0]*in.i[1];
		double v10=i[1]*in.i[0];
		if (v01<v10) return cInterval(v01,v11);
		return cInterval(v10,v11);
	}
	double v01=i[0]*in.i[1];
	double v10=i[1]*in.i[0];
	if (v01<v10) return cInterval(v01,v00);
	return cInterval(v10,v00);
	
}

cInterval cInterval::operator*(const double d) { 
	
	if (d>=0) return cInterval(i[0]*d,i[1]*d);
	else return  cInterval(i[1]*d,i[0]*d);

}

cInterval cInterval::operator/(const cInterval &in) { return *this*cInterval(1.0/in.i[1],1.0/in.i[0]); } // assumes that in does not contain 0

bool cInterval::operator^(const cInterval &in) { 
	
	// determines whether the intersection between the interval and 'in' is empty, 
	// returns true when it is, and compute it when it isn't
	
	if (i[1]<in.i[0]) return true;
	if (i[0]>in.i[1]) return true;
	if (i[1]>in.i[1]) i[1]=in.i[1];
	if (i[0]<in.i[0]) i[0]=in.i[0];
	return false;
	
}

cInterval cInterval::operator-() { return cInterval(-i[1],-i[0]); }

bool cInterval::isEmpty(double l, double r) { 
	
	// determines whether the intersection between the interval and [l,r] is empty, 
	
	if (i[1]<l) return true;
	if (i[0]>r) return true;
	return false;
	
}

double cInterval::getAbsLower() { 
	if (i[0]>=0) return i[0];
	if (i[1]>=0) return 0;
	return -i[1];
}

double cInterval::getAbsUpper() {
	if (i[0]+i[1]>=0) return i[1];
	return -i[0];
}

void cInterval::print() {

	std::cout << " [" << i[0] << "," << i[1] << "]" << std::endl;
	std::cout << std::endl;

}

cInterval operator*(double d, cInterval in) {

	if (d>=0) return cInterval(in.i[0]*d,in.i[1]*d);
	else return  cInterval(in.i[1]*d,in.i[0]*d);

}

cInterval operator+(double d, cInterval in) {

	return cInterval(d+in.i[0],d+in.i[1]);

}

bool cInterval::contains(double v) {

	if (v<i[0]) return false;
	if (v>i[1]) return false;
	return true;

}

void cInterval::bound(double v) {

	if (v<i[0]) i[0]=v;
	if (v>i[1]) i[1]=v;

}

void cInterval::bound(const cInterval &j) {

	if (j.i[0]<i[0]) i[0]=j.i[0];
	if (j.i[1]>i[1]) i[1]=j.i[1];

}

void cInterval::expand(double r) {

	i[0]-=r;
	i[1]+=r;

}

void cInterval::shrink(double r) {

	double _center=center();
	
	i[0]+=r;
	if (i[0] > _center) 
		i[0]=i[1]=_center;
	else 
		i[1] -= r;
	
}

