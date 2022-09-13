/*************************************************************************\
 
 NANO-D, 2010
 All Rights Reserved.
 
 \**************************************************************************/

#pragma once

#include <cstdlib>
#include <vector>
#include <map>

typedef struct {

	double radius;
	double mass;

} elementInformation;

class cPeriodicTable {
 public:
	static const std::map<int, elementInformation> periodicTableMap;
  // returns negative value for unknown elements
  static float getMass(const char element[2]);
  // returns negative value for unknown elements
  static float getRadius(const char element[2]);
};
