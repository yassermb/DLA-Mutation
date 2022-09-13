/*************************************************************************\

 Copyright 2019 Stephane Redon, Sergei Grudinin, CNRS / INRIA.
 All Rights Reserved.

 \**************************************************************************/


#include "cSpatialTransform.hpp"

#include "cVector3.hpp"

cSpatialTransform::cSpatialTransform()    { position.setZero();orientation.setIdentity(); }


cSpatialTransform::cSpatialTransform(cMatrix33 orient, cVector3 pos) : orientation(orient), position(pos) {}

void  cSpatialTransform::computeInverse(cSpatialTransform &inverse) const {

    // This function computes the inverse of the transform
    // A' = R*A + T
    // A = R^-1 A' - R^-1 T

    inverse.orientation = orientation.getTranspose();
    inverse.position = - (inverse.orientation*position);

    //    inverse.position.v[0]=-(inverse.orientation.m[0][0]*position.v[0]+inverse.orientation.m[0][1]*position.v[1]+inverse.orientation.m[0][2]*position.v[2]);
    //    inverse.position.v[1]=-(inverse.orientation.m[1][0]*position.v[0]+inverse.orientation.m[1][1]*position.v[1]+inverse.orientation.m[1][2]*position.v[2]);
    //    inverse.position.v[2]=-(inverse.orientation.m[2][0]*position.v[0]+inverse.orientation.m[2][1]*position.v[1]+inverse.orientation.m[2][2]*position.v[2]);

}

cSpatialTransform  cSpatialTransform::getInverse() const {

    cSpatialTransform inverse;
    // This function computes the inverse of the transform
    // A' = R*A + T
    // A = R^-1 A' - R^-1 T

    inverse.orientation = orientation.getTranspose();
    inverse.position = - (inverse.orientation*position);
    return inverse;
}

void cSpatialTransform::transformPoint(cVector3 &v, cVector3 &result) {
    result=orientation * v + position;
}

void cSpatialTransform::print()    {

    for (int i=0;i<3;i++) {

//        for (int j=0;j<3;j++) std::cout << "\t[" << orientation[3*i+j] << "]";
        for (int j=0;j<3;j++) std::cout << "\t[" << orientation.m[i][j] << "]";
        std::cout << "\t[" << position[i] << "]";
        std::cout << std::endl;

    }

    std::cout << std::endl;
    std::cout << std::endl;

}

