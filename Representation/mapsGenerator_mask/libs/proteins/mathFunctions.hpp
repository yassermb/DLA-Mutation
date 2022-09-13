#pragma once

#include "cVector3.hpp"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace math_functions {

// Compute angle ABC
// If angle is indeterminated, function returns false and doesn't change |angle|
// If successful, the variable |angle| gets the value of angle ABC in degrees.
bool computeAngle(const cVector3 &A, const cVector3 &B, const cVector3 &C, double *angle);

template <typename T>
int sgn(T val) { return (T(0) < val) - (val < T(0)); }

// Compute dihedral angle ABCD
// If angle is indeterminated, function returns false and doesn't change |dihAngle|
// If successful, |dihAngle| gets the value of angle ABCD in degrees [-180, 180].
bool computeDihedralAngle(const cVector3 &A, const cVector3 &B,
                          const cVector3 &C, const cVector3 &D, double *dihAngle);

void setDihedral(const cVector3 &A, const cVector3 &B, const cVector3 &C,
                 cVector3 *D, double CD, double B_C_D, double A_B_C_D);

} // math_functions
