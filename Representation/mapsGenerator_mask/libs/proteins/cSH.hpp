#pragma once

#include "cVector3.hpp"

typedef double Real;
typedef struct complex {
    Real x, y;
} sComplex;








class cSH {

    //constants
    const Real fourPi = 4*3.14159265358979323846;

    //math arrays
    Real     *factorial;
    Real     **Fcoeff;
    Real    **LegPoly;
    sComplex    **Y_C;
    sComplex    *Yxy;

    //aux init
    void makeFactorial(int P);
    void makeSphCoeff(int P);

    //aux functions
    void computeLegendre(const Real xval);
    void computeFourier(const Real  b);



    //memory function
    void allocHalfArray( Real** &A, int P);
    void dellocHalfArray( Real** &A, int P);
    void allocComplexArray( sComplex** &M, int P);
    void dellocComplexArray(sComplex** &M, int P);

    //parameters of the class
    int order;
    
public:
    cSH();
    ~cSH();
    sComplex    ** const getSH() const {return Y_C;}

    void init(int _P);
    static cVector3 cart2Sph( const cVector3 &v);
    void computeSpharmonic( const cVector3  &sv );


};
