#pragma once
#include "NeighbourSearch.hpp"
#include "cProtein.hpp"
#include <unordered_map>
#include <cMatrix33.hpp>
#include "cGrid.hpp"

// CSR and chosen residues
/////////////////////////////////////////////
typedef std::pair<size_t, char> _pair;
/////////////////////////////////////////////

class cProteinMapper {
protected:
		cProtein * _protein;
		uint8_t* _map;
		int32_t* _meta;
		const size_t gridSize;
		const float voxWidth;
		int* header;
public:
	int numRetype = -1;
	cProteinMapper(size_t gs, float vw ): gridSize(gs), voxWidth(vw) {}
	virtual ~cProteinMapper() {}
	virtual int getNbMaps(bool native, std::string scoreFileName) { return -1; }
	virtual int getMetaSize() { return -1; }
	virtual int getMapSize() { return -1; }
	virtual void setProtein(cProtein *protein) {
		_protein = protein; 
	};
	virtual void runCurrent() = 0;
	virtual void runGradient(float* grads, float* retypeMat, std::unordered_map<size_t, cVector3> &gradByAtom) = 0;
	virtual void next()=0;
	virtual const uint8_t* getMap();
	virtual int32_t* getMeta();
	virtual int* getHeader() { return header; };

	// CSR and chosen residues
	/////////////////////////////////////////////
	std::map<_pair, char> _mapSCR;
	std::map<_pair, char> _resList;
	bool _skipRes = false;
	cVector3 centerSphere;
	/////////////////////////////////////////////
};


class cResidueMapper : public cProteinMapper {
protected: 
	bool _orient;
	cVector3 origin;
	cVector3 vX;
	cVector3 vY;
	cVector3 vZ;
	cMatrix33 localToWorld;
	cMatrix33 worldToLocal;
	int currentChain;
	HighLevelIteratorGenerator<cProteinChain, cResidue>::iterator currentResIt;
	const float resSpread;
	virtual bool setLocalFrame(cResidue* res, cAtom * frameAtoms = NULL, bool orient = true);
    virtual bool setLocalFrameCa(cResidue* res, cAtom * frameAtoms = NULL, bool orient = true);
	virtual void setProtein(cProtein *protein);
	cResidueMapper(size_t gridSize, float resS, float voxWidth, bool orient = true) : cProteinMapper(gridSize, voxWidth),resSpread(resS), _orient(orient) {};
	virtual void next();
	virtual int getNbMaps(bool native, std::string scoreFileName);
};

class cAtomicResidueMapper : public cResidueMapper {
	int addAtom(cAtom*);
	cVector3 addGradAtom(cAtom* at , float* grads, float* typeVect);
	const float atomSpread;
	int nbType;
	const bool _hydrogens;
	const int _skipNeighb;
	static int armMetaSize(){return  16 * sizeof(int32_t);}

	int initResidueAtoms();
	int initResidueAtoms167();
	int initResidueAtoms4();
	int initResidueAtoms20();
	int initResidueAtoms21();
	std::map<cAminoResidue::eAminoType, std::map<std::string, int>> atomTypes;
	int getType(cAtom* atom);
	void runCurrentThread(int threadId);
	cGrid neighbourGrid;

public:
	cAtomicResidueMapper(size_t nt, size_t gridSize, float atomSpread, float voxWidth, bool hydrogen = false, bool orient = true, int skipNeighb = 0);
	virtual ~cAtomicResidueMapper();
	virtual void runCurrent();
	virtual void setProtein(cProtein *protein);
	virtual void runGradient(float* grads, float* retypeMat, std::unordered_map<size_t, cVector3> &gradByAtom);
	//virtual void runGradient2(float* grads, float* retypeMat, std::unordered_map<size_t, cVector3> &gradByAtom);
	//virtual const uint8_t* getMap();
	//virtual int32_t* getMeta();
	virtual int getMetaSize() { return armMetaSize(); }
	virtual int getMapSize() { return gridSize*gridSize*gridSize*(nbType+6); }   //+4 for SCR, +2 for R/L
	//virtual void setProtein(cProtein *protein);

	// CSR and chosen residues
	//////////////////////////////////////////
	void loadSCRfile();
	void loadResidueList();
	//////////////////////////////////////////
};


class cCoarseResidueMapper : public cResidueMapper {
	//bool setLocalFrame(cResidue* res, cAtom * frameAtoms = NULL, bool orient = true);
	//int addResidue(cResidue*);

	std::map<cResidue *, std::vector<cResidue *>> neighbourRes;
	const float residueSpread;
	int nbType;
	static int crmMetaSize() { return  16 * sizeof(int32_t); }

	int initResidueTypes();
	std::vector<int> residueTypes;
	int getType(cResidue* res);
	int addRes(cResidue*);

public:
	cCoarseResidueMapper(size_t gridSize, float atomSpread, float voxWidth,  bool orient = true);
	virtual ~cCoarseResidueMapper();
	virtual void runCurrent();
	virtual void setProtein(cProtein *protein);
	virtual void runGradient(float* grads, float* retypeMat, std::unordered_map<size_t, cVector3> &gradByAtom) {};
	//virtual void runGradient2(float* grads, float* retypeMat, std::unordered_map<size_t, cVector3> &gradByAtom) {};
	//virtual const uint8_t* getMap();
	//virtual int32_t* getMeta();
	//virtual int getNbMaps(bool native, std::string scoreFileName);
	virtual int getMetaSize() { return crmMetaSize(); }
	virtual int getMapSize() { return gridSize*gridSize*gridSize*nbType*6; }
	//virtual void setProtein(cProtein *protein);
};
