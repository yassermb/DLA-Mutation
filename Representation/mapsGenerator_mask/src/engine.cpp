/*************************************************************************\

Sergei Grudinin, 2012
Mikhail Karasikov, 2016
All Rights Reserved.

\**************************************************************************/

#include "engine.hpp"

//#include "bessel.h" // for bessel functions
//#include "besselfrac.h" // for bessel functions
#include <assert.h>
#include <string.h>
#include <cstdlib>
#include <cstdio>

#include "cClassicalParserPDB.hpp"
#include "cGrid.hpp"
#include "cAlgorithmTemplates.hpp"
#include "cRotamerLibrary.hpp"
#include "mathFunctions.hpp"
#include "cSparseMatrixOutput.hpp"
#include "energyModel.hpp"
#include "cProteinFeaturizer.hpp"
#include "cProteinMapper.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <unordered_map>
#include "npy.hpp"

#include "cSH.hpp" //for SH and oriented graphs

const size_t kMaxRotamersPerResidue = 5;
const size_t kNumIterations = 10;

using std::vector;
using std::set;


bool write_vector(const Eigen::SparseVector<double> &vector,
                  const std::string &filename) {
  cSparseMatrixOutputMAT matrixOutput(vector.size(), 1);
  if (!matrixOutput.initialize(filename))
    return false;

  for (Eigen::SparseVector<double>::InnerIterator it(vector); it; ++it) {
    if (!matrixOutput.writeTriplet(it.index(), 0, it.value()))
      return false;
  }
  return matrixOutput.deinitialize();
}


bool writeMap(std::fstream & binOutputFile, const uint8_t* map, const int32_t* meta, size_t mapSize, size_t metaSize) {
	binOutputFile.write((char*)map, mapSize);
	binOutputFile.write((char*)meta, metaSize);
	return true;
}

bool writeMap(FILE* binOutputFile, const uint8_t* map, const int32_t* meta, size_t mapSize, size_t metaSize) {
	fwrite(map, sizeof(char), mapSize, binOutputFile);
	fwrite(meta, sizeof(char), metaSize, binOutputFile);
	return true;
}

bool writeSHfeatures(cProtein *protein, const std::string &outputFilename, const int order) {

    //allocate memory + precompute constants
    cSH sh;
    sh.init(order);

    //first, compute all local frames
    for (auto &residue : protein->residues()) {
        if(residue.setLocalFrame()==0) {
            std::cerr << "No backbone atoms in " << residue.resName << " "<< residue.seqNumber <<"\n";
            return false;
        }
    }

    //now, compute all vs all orientations
    printf("#printing SH for %lu residues and %d expansion order. In total %d x %lu values\n", protein->numResidues(), order, order*order, protein->numResidues()*(protein->numResidues()-1));
    for (auto &resi : protein->residues()) {
        for (auto &resj : protein->residues()) {

            if(&resi == &resj) { //diagonal case
                continue;
            }

            cVector3 worldJ = resj.getAtom("CA")->getPosition();
            cVector3 jInI = resi.worldToLocal*worldJ;

            cVector3 localSphCoords = cSH::cart2Sph(jInI);
            sh.computeSpharmonic(localSphCoords);

            //OK, now, SH coeeficients are in Y_C[n][m].x and Y_C[n][m].y!
            //We have order^2 values
            for (int l=0; l < order; l++) {
                for (int m=1; m <= l; m++) {
                    printf("%8.4f", -2*sh.getSH()[l][m].y); // the -m part
                }

                printf("%8.4f", sh.getSH()[l][0].x);

                for (int m=1; m <= l; m++) {
                    printf("%8.4f", +2*sh.getSH()[l][m].x); // the +m part
                }
            }
            printf("\n");
        }
    }

    return true;
}

bool writeAllMaps(cProtein *protein, cProteinMapper *mapper, const std::string &outputFilename, bool native, const std::string &scoreFilename, std::map<int, char> ssByRes, std::map<int, float> areaByRes) {
	std::fstream binOutputFile;
	mapper->setProtein(protein);
	binOutputFile.open(outputFilename, std::ios_base::out | std::ios_base::binary);
	FILE* binFile = fopen(outputFilename.c_str(),"wb");
	int nbMaps = mapper->getNbMaps(native, scoreFilename);
	int mapSize = mapper->getMapSize();
	int metaSize = mapper->getMetaSize();
	int* header = mapper->getHeader();
	header[7] = nbMaps;

  std::cout<<"nbMaps\t"<<nbMaps<<std::endl;

	//binOutputFile.write((char*) header, header[0]);
	fwrite(header, sizeof(char), header[0], binFile);
	std::ifstream scoreFile;
	std::string line;
	if (!scoreFilename.empty()) {
		scoreFile.open(scoreFilename, std::ios_base::in);
	}

	int resCode = -1;
	float score = -1;
	bool eof = false;
	int nbMappedRes = 0;
	
		for (int indexMap = 0; indexMap < protein->numResidues(); indexMap++) {
			mapper->runCurrent();
			const uint8_t* map = mapper->getMap();
			int32_t* meta = mapper->getMeta();
			if (!scoreFilename.empty()) {
				while (resCode < ((int*)meta)[0] && !eof && scoreFile) {
					std::getline(scoreFile, line);
					if (line.empty()) {
						eof = true;
						break;
					}
					int count = 0;
					while (line.at(count) != 'r') {
						count++;
					}
					assert(line.at(count) == 'r');
					count++;
					assert(line.at(count) == '<');
					count++;
					int intStart = count;
					while (line.length() > count && line.at(count) != '>') {
						count++;
					}
					resCode = stoi(line.substr(intStart, count - intStart));
					assert(line.at(count) == '>');
					count++;
					if (line.at(count) != 'R'){
                                                count+=4;
                                        }

					assert(line.at(count) == 'R');
					count++;
					assert(line.at(count) == '<');
					while (line.length() > count && line.at(count) != '>') {
                                                count++;
                                        }
					assert(line.at(count) == '>');
					count++;
					assert(line.at(count) == ' ');
					count++;
					score = stof(line.substr(count));
				}
			}
			int scoreInt;
			if (native) {
				scoreInt = 1000000;
			}
			else if (resCode == ((int*)meta)[0]) {
				scoreInt = 1000000 * score;
			}
			else {
				scoreInt = -1000000;
			}
			meta[2] = scoreInt;
			if (ssByRes.find(((int*)meta)[0]) != ssByRes.end()) {
				meta[3] = ssByRes.at(((int*)meta)[0]);
				meta[4] = 10 * areaByRes.at(((int*)meta)[0]);
			}
			if (scoreInt >= 0 && !mapper->_skipRes) {
				//bool wm = writeMap(binOutputFile, map, meta, mapSize, metaSize);
				bool wm = writeMap(binFile, map, meta, mapSize, metaSize);
				
				if (!wm) return wm;
			}
			mapper->next();
			nbMappedRes++;
		}
		
		//std::cout << "Nb mapped residues : " << nbMappedRes << std::endl;
	return true;
}


bool computeGrad(cProtein *protein, cProteinMapper *mapper, const std::string &gradientMap, const std::string &outputFilename) {
	std::unordered_map<size_t, cVector3> gradByAtom;
	for (cAtom at : protein->atoms()) {
		gradByAtom.insert({ {at.atomSerial, cVector3(0, 0, 0)} });
	}
	protein->complete(true);
	mapper->setProtein(protein);

	int nbMaps = mapper->getNbMaps(true, "");
	int mapSize = mapper->getMapSize();
	std::vector<unsigned long> retypeMatShape;
	std::vector<float> retypeMat;

	std::ifstream stream(gradientMap, std::ifstream::binary);
	if (!stream) {
		throw std::runtime_error("io error: failed to open a file.");
	}
	{
		std::string header = npy::read_header(stream);
		// parse header
		bool fortran_order;
		std::string typestr;
		npy::parse_header(header, typestr, fortran_order, retypeMatShape);

		// check if the typestring matches the given one
		npy::Typestring typestring_o{ retypeMat };
		std::string expect_typestr = typestring_o.str();
		if (typestr != expect_typestr) {
			throw std::runtime_error("formatting error: typestrings not matching");
		}

		// compute the data size based on the shape
		auto size = static_cast<size_t>(npy::comp_size(retypeMatShape));
		retypeMat.resize(size);

		// read the data
		stream.read(reinterpret_cast<char*>(retypeMat.data()), sizeof(float)*size);
	}

	mapper->numRetype = retypeMatShape[1];

	int resCode = -1;
	float score = -1;
	bool eof = false;
	for (int indexMap = 0; indexMap < protein->numResidues(); indexMap++) {
		std::vector<float> gradArray;
		std::vector<unsigned long> gradShape;
		{
			std::string header = npy::read_header(stream);
			// parse header
			bool fortran_order;
			std::string typestr;
			npy::parse_header(header, typestr, fortran_order, gradShape);

			// check if the typestring matches the given one
			npy::Typestring typestring_o{ gradArray };
			std::string expect_typestr = typestring_o.str();
			if (typestr != expect_typestr) {
				throw std::runtime_error("formatting error: typestrings not matching");
			}

			// compute the data size based on the shape
			auto size = static_cast<size_t>(npy::comp_size(gradShape));
			gradArray.resize(size);

			// read the data
			stream.read(reinterpret_cast<char*>(gradArray.data()), sizeof(float)*size);
		}
		mapper->runGradient(gradArray.data(),retypeMat.data(), gradByAtom);
		
		mapper->next();
	}


	cVector3 sumGradAll(0, 0, 0);


	FILE* f = fopen(outputFilename.c_str(), "w");
	int n = 0;
	for (cAtom at : protein->atoms()) {
		n++;
		if (gradByAtom.find(at.atomSerial) != gradByAtom.end()) {
			at.setPosition(gradByAtom[at.atomSerial]);
			sumGradAll += gradByAtom[at.atomSerial];
			gradByAtom.erase(at.atomSerial);
			at.toFile(f);
		}
		else {

			at.setPosition(cVector3(0,0,0));
			at.toFile(f);
		}
	}

	return true;
}


bool writeFeatureVector(cProtein *protein,
                        cProteinFeaturizer *featurizer,
                        const std::string &filename) {
  auto feature_vector = featurizer->featurize(protein);
  return write_vector(feature_vector, filename);
}

void greedySearch(cProtein *protein, size_t numIterations) {
  printf("Iteration     Energy           Changed rotamers\n");

  vector<size_t> rotamers;
  for (auto &residue : protein->residues()) {
    cAminoResidue *res = dynamic_cast<cAminoResidue *>(&residue);
    if (!res || res->numRotamers() < 2)
      continue;

    res->setRotamer(0);
    rotamers.push_back(0);
  }

  energy::initializeNeighbourAtoms(protein);
  for (size_t k = 0; k < numIterations; ++k) {

    double totalEnergySum = 0;
    size_t numChangedRotamers = 0;
    size_t idx = 0;
    for (auto &residue : protein->residues()) {
      cAminoResidue *res = dynamic_cast<cAminoResidue *>(&residue);
      if (!res || res->numRotamers() < 2)
        continue;

      double minEnergy = 1E20;
      size_t bestRotamer = -1;
      double libEnergy;

      for (size_t j = 0; j < res->numRotamers() && j < kMaxRotamersPerResidue; ++j) {
        res->setRotamer(j, &libEnergy);
        double energy = energy::totalEnergy(*res);
        energy += libEnergy;
        if (energy < minEnergy) {
          minEnergy = energy;
          bestRotamer = j;
        }
      }
      assert(bestRotamer != static_cast<size_t>(-1));
      if (bestRotamer != rotamers[idx]) {
        numChangedRotamers++;
        rotamers[idx] = bestRotamer;
      }
      res->setRotamer(bestRotamer);
      totalEnergySum += minEnergy;
      idx++;
    }
    printf("%9zd     %6lf          %zd\n", k, totalEnergySum, numChangedRotamers);
    if (!numChangedRotamers)
      break;
  }
#if 0
  size_t idx = 0;
  for (residueIterator resIt = protein->residuesBegin();
                        resIt != protein->residuesEnd(); ++resIt) {
    cAminoResidue *res = dynamic_cast<cAminoResidue *>(&(*resIt));
    if (!res || res->numRotamers() < 2)
      continue;
    for (size_t i = 0; i < res->numRotamers() && i < kMaxRotamersPerResidue; ++i) {
      if (i == rotamers[idx]) {
        printf("1 ");
      } else {
        printf("0 ");
      }
    }
    idx++;
  }
  printf("\n%d\n", rotamers.size());
#endif
}

void doSmth(cProtein *protein) {
  // writeEnergyMatrix("Energy");
  // setRotamers("vector.txt");

	// some output
#if 1
  greedySearch(protein, kNumIterations);
#else
  energy::initializeNeighbourAtoms(protein);

  for (size_t k = 1; k <= kNumIterations; ++k) {

    for (residueIterator resIt = protein->residuesBegin();
                          resIt != protein->residuesEnd(); ++resIt) {
      cAminoResidue *res = dynamic_cast<cAminoResidue *>(&(*resIt));
      if (!res || !res->numRotamers())
        continue;

      double min_energy = 1000000;
      double initial_energy = energy::totalEnergy(*resIt);
      double libEnergy;
      size_t best_rotamer;
      for (size_t j = 0; j < res->numRotamers(); ++j) {
        res->setRotamer(j, &libEnergy);
        double energy = energy::totalEnergy(*res);
        energy += libEnergy;
        if (energy < min_energy) {
          min_energy = energy;
          best_rotamer = j;
        }
      }
      double prob = res->setRotamer(best_rotamer, &libEnergy);
      if (k == kNumIterations) {
        std::cout << "Chain " << res->chainId << ": Residue "
                  << res->seqNumber << ": phi " << res->getPhi()
                  << ", psi " << res->getPsi()
                  << " is set to " << best_rotamer
                  << " rotamer with probability: " << prob
                  << " Total: " << min_energy
                  << " Initial: " << initial_energy
                  << " Relative: " << libEnergy << std::endl;
      }
    }
  }
#endif
/*
  cChain::atomIterator atom = chains[0]->atomsBegin();
  for (set<cAtom *>::iterator it = atom->neighbours.begin();
                                     it != atom->neighbours.end(); ++it)
    std::cout << *it << std::endl;
*/
}

bool writeEnergyMatrix(cProtein *protein, const std::string &filename) {
  vector<cAminoResidue *> residues;
  vector<vector<cAminoResidue *>> neighbours
      = energy::initializeNeighbourResidues(protein, &residues);
  // reindexing residues and counting number of all rotamers
  size_t numRotamers = 0;
  for (size_t resIdx = 0; resIdx < residues.size(); ++resIdx) {
    residues[resIdx]->residueIndex = numRotamers;
    numRotamers += std::min(residues[resIdx]->numRotamers(),
                            kMaxRotamersPerResidue);
  }
#if 0
  cSparseMatrixOutputASCII matrixOutput;
#else
  cSparseMatrixOutputMAT matrixOutput(numRotamers, numRotamers + 3);
#endif
  if (!matrixOutput.initialize(filename))
    return false;

  // Write inexes of all rotamers
  // residues
  for (size_t resIdx = 0; resIdx < residues.size(); ++resIdx) {
    cAminoResidue &residue = *residues[resIdx];

    // rotamers
    for (size_t rotIdx = 0; rotIdx < kMaxRotamersPerResidue
                            && rotIdx < residue.numRotamers(); ++rotIdx) {

      if (!matrixOutput.writeTriplet(residue.residueIndex + rotIdx, 0,
                                      residue.residueIndex))
        return false;
    }
  }
  // Write libEnergy
  // residues
  for (size_t resIdx = 0; resIdx < residues.size(); ++resIdx) {
    cAminoResidue &residue = *residues[resIdx];

    // rotamers
    for (size_t rotIdx = 0; rotIdx < kMaxRotamersPerResidue
                            && rotIdx < residue.numRotamers(); ++rotIdx) {

      double libEnergy;
      residue.setRotamer(rotIdx, &libEnergy);
      if (!matrixOutput.writeTriplet(residue.residueIndex + rotIdx, 1,
                                      libEnergy))
        return false;
    }
  }
  // Write frameEnergy
  // residues
  for (size_t resIdx = 0; resIdx < residues.size(); ++resIdx) {
    cAminoResidue &residue = *residues[resIdx];

    // rotamers
    for (size_t rotIdx = 0; rotIdx < kMaxRotamersPerResidue
                            && rotIdx < residue.numRotamers(); ++rotIdx) {
      residue.setRotamer(rotIdx);
      if (!matrixOutput.writeTriplet(residue.residueIndex + rotIdx, 2,
                                      energy::frameEnergy(residue)))
        return false;
    }
  }
  // Write pairEnergy
  // residues
  for (size_t resIdx = 0; resIdx < residues.size(); ++resIdx) {
    cAminoResidue &residue = *residues[resIdx];

    // rotamers
    for (size_t rotIdx_1 = 0; rotIdx_1 < kMaxRotamersPerResidue
                              && rotIdx_1 < residue.numRotamers(); ++rotIdx_1) {
      residue.setRotamer(rotIdx_1);
      if (!matrixOutput.writeTriplet(residue.residueIndex + rotIdx_1,
                                      3 + residue.residueIndex + rotIdx_1,
                                      energy::pairEnergy(residue, residue) / 2))
        return false;
      // neighbour residues
      for (size_t nbIdx = 0; nbIdx < neighbours[resIdx].size(); ++nbIdx) {
        cAminoResidue &neighbour = *neighbours[resIdx][nbIdx];

        // rotamers of neighbour residue
        for (size_t rotIdx_2 = 0; rotIdx_2 < kMaxRotamersPerResidue
                                  && rotIdx_2 < neighbour.numRotamers(); ++rotIdx_2) {

          neighbour.setRotamer(rotIdx_2);
          if (!matrixOutput.writeTriplet(neighbour.residueIndex + rotIdx_2,
                                          3 + residue.residueIndex + rotIdx_1,
                                          energy::pairEnergy(residue, neighbour)))
            return false;
        }
      }
    }
  }
  return matrixOutput.deinitialize();
}

inline bool areAnglesEqual(double angle_1, double angle_2, double dev) {
  double diff = std::abs(angle_1 - angle_2);
  return std::abs(diff) < dev || std::abs(diff - 360.0) < dev;
}

double chi1Quality(const cProtein &protein, const cProtein &ethalon) {
  size_t numAngles = 0;
  size_t numEqualAngles = 0;

  for (auto resIt = protein.residues().begin(), ethalonResIt = ethalon.residues().begin();
            resIt != protein.residues().end() && ethalonResIt != ethalon.residues().end();
                ++resIt, ++ethalonResIt) {
    const cAminoResidue *res
      = dynamic_cast<const cAminoResidue *>(&(*resIt));
    const cAminoResidue *ethalonRes
      = dynamic_cast<const cAminoResidue *>(&(*ethalonResIt));
    if (!res)
      continue;
    numAngles++;
    numEqualAngles += areAnglesEqual(res->getChi(1), ethalonRes->getChi(1), 40);
  }
  return static_cast<double>(numEqualAngles) / numAngles * 100;
}

double chi12Quality(const cProtein &protein, const cProtein &ethalon) {
  size_t numAngles = 0;
  size_t numEqualAngles = 0;

  for (auto resIt = protein.residues().begin(), ethalonResIt = ethalon.residues().begin();
            resIt != protein.residues().end() && ethalonResIt != ethalon.residues().end();
                ++resIt, ++ethalonResIt) {
    const cAminoResidue *res
      = dynamic_cast<const cAminoResidue *>(&(*resIt));
    const cAminoResidue *ethalonRes
      = dynamic_cast<const cAminoResidue *>(&(*ethalonResIt));
    if (!res)
      continue;
    numAngles++;
    numEqualAngles += areAnglesEqual(res->getChi(1), ethalonRes->getChi(1), 40)
                      * areAnglesEqual(res->getChi(2), ethalonRes->getChi(2), 40);
  }
  return static_cast<double>(numEqualAngles) / numAngles * 100;
}

double rmsdQuality(const cProtein &protein, const cProtein &ethalon) {
  double squareSum = 0.0;
  size_t numAtoms = 0;

  for (auto resIt = protein.residues().begin(),
            ethalonResIt = ethalon.residues().begin();
                resIt != protein.residues().end() &&
                ethalonResIt != ethalon.residues().end();
                    ++resIt, ++ethalonResIt) {
    for (const auto &first_atom : resIt->atoms()) {
      for (const auto &second_atom : ethalonResIt->atoms()) {
        if (!strcmp(first_atom.name, second_atom.name)) {
          numAtoms++;
          squareSum += (first_atom.getPosition() - second_atom.getPosition()).norm2();
        }
      }
    }
  }
  return sqrt(squareSum / numAtoms);
}

size_t clashQuality(const cProtein &protein) {
  return energy::numClashes(protein);
}

double energyQuality(cProtein *protein) {
  double total_energy = 0;

  energy::initializeNeighbourAtoms(protein);

  for (const auto &residue : protein->residues()) {
    if (const auto *amino_residue = dynamic_cast<const cAminoResidue *>(&residue))
      total_energy += energy::totalEnergy(*amino_residue);
  }
  return total_energy / 2;
}
