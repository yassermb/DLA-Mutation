/*************************************************************************\

 Copyright 2019 Sergei Grudinin, Stephane Redon, CNRS / INRIA.
 All Rights Reserved.

 \**************************************************************************/


#pragma once


#include "cVector3.hpp"
#include "cMatrix33.hpp"

#define CSPATIALTRANSFORM_PI                            3.1415926535897932
#define CSPATIALTRANSFORM_HALF_PI                        1.57079632679489661

// ***************************************************************************************
// ***************************************************************************************
// **                                                                                    **
// ** class : cSpatialTransform                                                            **
// **                                                                                    **
// ** Spatial Transform information                                                        **
// **                                                                                    **
// ***************************************************************************************
// ***************************************************************************************

// This class handles reference frames, and includes some useful functions to treat the transform
// as a regular 3d transform or as a spatial transform.

class cSpatialTransform {

    public :

    cMatrix33                orientation;
    cVector3                position;

    cSpatialTransform();
    cSpatialTransform(cMatrix33 orient, cVector3 pos);

    void                    computeInverse(cSpatialTransform &inverse) const;
    cSpatialTransform       getInverse() const;
    void                    rightMultiply3d(cSpatialTransform &leftMember, cSpatialTransform &result);
    void                    rightMultiplyInverse3d(cSpatialTransform &leftMember, cSpatialTransform &result);
    void                    transformPoint(cVector3 &v, cVector3 &result);
    cVector3                operator*(const cVector3 &v) const { return orientation*v+position; }
    cSpatialTransform        operator*(const cSpatialTransform &t) const { return cSpatialTransform(orientation*t.orientation,orientation*t.position+position); }
    cMatrix33                operator*(const cMatrix33 &t) const { return orientation*t; }

    void                    print();

    double                    test() {cVector3 tmp(1.0,1.0,1.0); return (orientation*tmp).norm();}

    //    cIAVector3                operator*(const cIAVector3 &v) const;
    cSpatialTransform        operator-(cSpatialTransform &t)  {
        return cSpatialTransform(orientation - t.orientation,position - t.position);
    }

};

typedef cSpatialTransform *pcSpatialTransform;


inline std::ostream& operator<<(std::ostream& os, cSpatialTransform &vec) { os << "";vec.print(); return os;}

