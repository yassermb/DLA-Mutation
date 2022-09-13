#include "cSH.hpp"
#include <math.h>       /* sqrt */

cSH::cSH() {
    
    // math arrays
    factorial = NULL;
    Fcoeff = NULL;
    Yxy = NULL;
    LegPoly = NULL;
    Y_C = NULL;
}

cSH::~cSH() {

    // math arrays
    if (factorial) delete [] factorial;

    if (Fcoeff) dellocHalfArray(Fcoeff, order);

    if (Yxy) delete [] Yxy;

    if (LegPoly) dellocHalfArray(LegPoly, order);

    if (Y_C) dellocComplexArray(Y_C, order);

}

void cSH::init(int _P) {

    order = _P;

    factorial=new Real[2*order+2];
    makeFactorial(order);

    allocHalfArray(Fcoeff, order);
    makeSphCoeff(order);

    allocComplexArray(Y_C, order);

    Yxy=new sComplex [order];
    allocHalfArray(LegPoly, order);

}

//FIXME - at low rho Maclaurin series should be used instead!!!
cVector3 cSH::cart2Sph( const cVector3 &v)
{
    cVector3 s;

    /* rho */
    s[0]=v.norm();

    if (s[0] < fabs(v[2])) s[0]=fabs(v[2]);

    /* alpha (Theta) */
    if (s[0]==0.0) s[1]=0.0;
    else s[1]=acos(v[2]/s[0]);

    /* beta (Phi) */
    if ((v[0]==0.0) && (v[1]==0.0)) s[2]=0.0;
    else s[2]=atan2(v[1], v[0]);

    return s;
} /* Cart2Sph */

void cSH::computeSpharmonic( const cVector3  &sv ) {

    int     m, n;
    Real    ytemp;

    computeLegendre(cos(sv[1]));
    computeFourier(sv[2]);

    //order expansion order
    for (n=0; n < order; n++) {
        for (m=0; m <= n; m++) {
            ytemp=  Fcoeff[n][m] * LegPoly[n][m];
            Y_C[n][m].x=ytemp * Yxy[m].x;
            Y_C[n][m].y=ytemp * Yxy[m].y;
        } /* for m */
    } /* for n */
} /* computeSpharmonic */

void cSH::computeLegendre(const Real xval) {

    int Li, Lj;
    Real negterm, oddfact, nextodd, sqroot, sqrtterm;

    negterm=1.0;
    oddfact=1.0;
    nextodd=1.0;
    sqroot=sqrt(1.0 - xval*xval);
    sqrtterm=1.0;
    for (Li=0;Li < order;Li++) {
        LegPoly[Li][Li]=negterm*oddfact*sqrtterm;
        negterm=-negterm;
        oddfact *= nextodd;
        nextodd += 2.0;
        sqrtterm *= sqroot;
        if (Li < order-1) {
            LegPoly[Li+1][Li]=xval * (Real)(2*Li+1) * LegPoly[Li][Li];
            for (Lj=Li+2;Lj < order;Lj++) {
                LegPoly[Lj][Li]=(xval*(Real)(2*Lj-1)*LegPoly[Lj-1][Li] -
                                 (Real)(Lj+Li-1)*LegPoly[Lj-2][Li])/(Real)(Lj-Li);
            }
        }
    }
}

void cSH::computeFourier(const Real  b) {

    int   m;

    for (m=0; m < order; m++) {
        Yxy[m].x=cos((Real) m * b);
        Yxy[m].y=sin((Real) m * b);
    } /* for m */

} /* computeFourier */

void cSH::makeFactorial(int P) {

    int n;
    factorial[0]=1.0;
    for (n=1; n < 2 * (P + 1); n++) {
        factorial[n]=(Real) n *factorial[n - 1];
    } /* for n */
}

void cSH::makeSphCoeff(int P) {
    int l, m;
    for (l=0; l < P; l++) {
        Real coeff = (2*l+1)/fourPi;
        for (m=0; m <= l; m++) {
            Fcoeff[l][m]= sqrt(coeff*factorial[l-m]/factorial[l+m]);
        } /* for m */
    } /* for l */
}
void cSH::allocHalfArray( Real** &A, int P) {

    A=new Real* [P];
    for (int i=0;i<P;i++)
        *(A+i)=new Real[i+1];
}

void cSH::dellocHalfArray( Real** &A, int P) {

    for (int i=0;i<P;i++)
        delete [] A[i];
    delete [] A;
}

void cSH::allocComplexArray( sComplex** &M, int P) {

    M=(sComplex**) new sComplex* [P];
    for (int i=0;i<P;i++)
        *(M+i)=new sComplex[i+1];
}

void cSH::dellocComplexArray(sComplex** &M, int P) {

    for (int i=0;i<P;i++)
        delete [] M[i];
    delete [] M;
}

