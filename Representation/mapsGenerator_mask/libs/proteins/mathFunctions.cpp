#include "mathFunctions.hpp"

#include <assert.h>


bool math_functions::computeAngle(const cVector3 &A, const cVector3 &B,
                                  const cVector3 &C, double *angle) {

  if ((A - B).isZero() || (C - B).isZero())
    return false;

  *angle = acos((A - B).normalize() | (C - B).normalize()) * 180.0 / M_PI;
  return true;
}

bool math_functions::computeDihedralAngle(const cVector3 &A, const cVector3 &B,
                                          const cVector3 &C, const cVector3 &D,
                                          double *dihAngle) {
  cVector3 n1 = (B - A) ^ (C - B);
  cVector3 n2 = (C - B) ^ (D - C);
  if (n1.isZero() || n2.isZero())
    return false;

  *dihAngle = acos(n1.normalize() | n2.normalize()) * sgn(n1 | (D - C)) * 180.0 / M_PI;
  return true;
}

void math_functions::setDihedral(const cVector3 &A, const cVector3 &B, const cVector3 &C,
                                 cVector3 *D, double CD, double B_C_D, double A_B_C_D) {
  assert(D);

  cVector3 e_1 = C - B;
  e_1.normalize();
  cVector3 e_2 = A - B;
  e_1.orthogonalize(&e_2);
  e_1.rotate(&e_2, A_B_C_D * M_PI / 180.0);
  e_2.normalize();
  *D = C + (-e_1 * cos(B_C_D * M_PI / 180.0)
            + e_2 * sin(B_C_D * M_PI / 180.0)) * CD;
}
