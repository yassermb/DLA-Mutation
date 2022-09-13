/*************************************************************************\

Sergei Grudinin, 2012
Mikhail Karasikov, 2016
All Rights Reserved.

\**************************************************************************/


#pragma once
#include <string>

#include "cProtein.hpp"
#include "cProteinFeaturizer.hpp"
#include "cProteinMapper.hpp"


bool writeFeatureVector(cProtein *protein,
                        cProteinFeaturizer *featurizer,
                        const std::string &);

bool writeEnergyMatrix(cProtein *protein, const std::string &);

bool writeSHfeatures(cProtein *protein, const std::string &outputFilename, const int order);

bool writeAllMaps(cProtein *protein, cProteinMapper *mapper, const std::string &outputFilename, bool native, const std::string &scoreFilename , std::map<int, char> ssByRes, std::map<int, float> areaByRes);

bool computeGrad(cProtein *protein, cProteinMapper *mapper, const std::string &gradientMap, const std::string &outputFilename);

void greedySearch(cProtein *protein, size_t numIterations);

// Root mean square deviation
double rmsdQuality(const cProtein &protein, const cProtein &ethalon);
// Percent correct of chi_1
// (within 40 degrees those of the native)
double chi1Quality(const cProtein &protein, const cProtein &ethalon);
// Percent correct of chi_1 and chi_2
// (within 40 degrees those of the native)
double chi12Quality(const cProtein &protein, const cProtein &ethalon);
// Number of atom pairs in clash
size_t clashQuality(const cProtein &protein);
// Energy of the protein
double energyQuality(cProtein *protein);
// Some stuff
void doSmth(cProtein *protein);
