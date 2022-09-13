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

// The class cInterval implements the basic interval operations on real intervals


class cInterval {

public:

	double i[2];

	cInterval();
	cInterval(double v);
	cInterval(double ll, double rr);

	cInterval	operator+(const cInterval &in);
	cInterval	operator-(const cInterval &in);
	cInterval&	operator+=(const cInterval &in);
	cInterval&	operator+=(const double d);
	cInterval&	operator-=(const cInterval &in);
	cInterval	operator*(const cInterval &in);
	cInterval	operator*(const double d);
	cInterval	operator/(const cInterval &in);
	bool		operator^(const cInterval &in);
	cInterval	operator-();

	bool		isEmpty(double l, double r);
	double		getAbsLower();
	double		getAbsUpper();
	bool		contains(double v);
	void		bound(double v);
	void		bound(const cInterval &j);
	void		expand (double r);
	void		shrink (double r);

	void		print();
	double		center() { return 0.5*(i[0]+i[1]); }
	double		diameter() { return i[1]-i[0]; }

};

cInterval operator*(double d, cInterval in);
cInterval operator+(double d, cInterval in);


