#include "cProteinMapper.hpp"
#include "cAtom.hpp"
#include "cGrid.hpp"
#include "cAminoResidue.hpp"
#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <assert.h>
#include <thread>
#include <math.h>

// SCR and chosen residues 
/////////////////////////////////////////////
#include <string>
#include <sstream>
#include <map>
/////////////////////////////////////////////

#include<chrono>

cAtomicResidueMapper::cAtomicResidueMapper(size_t nt , size_t gs, float as, float vw, bool h, bool orient, int removeNeighb): cResidueMapper(gs, as, vw, orient), atomSpread(as), _hydrogens(h), nbType(nt), _skipNeighb(removeNeighb)
{

	nbType = initResidueAtoms();
	int metaSize = armMetaSize();
	_meta = new int32_t[metaSize/sizeof(int32_t)];
	memset(_meta, 0, metaSize );
	_map = new uint8_t[gs*gs*gs*(nbType+6)]; //+4 for SCR, +2 for R/L
	int headerSize = 8;
	header = new int[headerSize];
	header[0] = headerSize * sizeof(int);
	header[1] = 7919;
	header[2] = nbType+6; //+4 for SCR, +2 for R/L
	header[3] = gridSize;
	header[4] = gridSize;
	header[5] = gridSize;
	header[6] = metaSize;
	header[7] = -1;

	// SCR and chosen residues
	/////////////////////////////////////////////
	loadSCRfile();
	for(auto elem : _mapSCR)
	{
   		std::cout << elem.first.first << "\t" << elem.first.second << "\t" << elem.second << "\n";
	}
	loadResidueList();
	for(auto elem : _resList)
	{
   		std::cout << elem.first.first << "\t" << elem.first.second << "\t" << elem.second << "\n";
	}
	/////////////////////////////////////////////
}

// SCR
/////////////////////////////////////////////
void cAtomicResidueMapper::loadSCRfile(){
	std::string line;
	std::ifstream transform_file ("scrinfo");
	int resNum;
	char chID, scr;
	if (transform_file.is_open())
	{
		while(getline (transform_file,line))
		{    
			std::stringstream token_line (line); 
			std::string intermediate;
			getline (token_line, intermediate, ';');
			resNum = stoi(intermediate);
			getline (token_line, intermediate, ';');
			chID = intermediate[0];
			getline (token_line, intermediate, ';');
			scr = intermediate[0];
			// std::cout<<"intermediate:\t"<<intermediate<<std::endl;
			_mapSCR.insert({{resNum, chID}, scr});
		}
		transform_file.close();
	}
	else 
		printf("Unable to open file");
}
/////////////////////////////////////////////


// Chosen residues
/////////////////////////////////////////////
void cAtomicResidueMapper::loadResidueList(){
	std::string line;
	std::ifstream transform_file ("resinfo");
	int resNum;
	char chID;
	if (transform_file.is_open())
	{
		while(getline (transform_file,line))
		{    
			std::stringstream token_line (line); 
			std::string intermediate;
			getline (token_line, intermediate, ';');
			resNum = stoi(intermediate);
			getline (token_line, intermediate, ';');
			chID = intermediate[0];
			_resList.insert({{resNum, chID}, chID});
		}
		transform_file.close();
	}
	else 
		printf("Unable to open file");
}
/////////////////////////////////////////////


#define DEFINE_ATOMS(...) { static std::string cl[] = {__VA_ARGS__}; \
                            for (size_t i = 0; i < sizeof(cl) / sizeof(std::string); ++i){ \
								name_type_map.insert(std::pair<std::string, int>(cl[i], n));\
								/*std::cout << cr << " "<< cl[i]<< std::endl;*/\
							n++;}}
#define DEFINE_ALIASE(al_1, al_2) name_type_map.insert(std::pair<std::string, int>(al_1, name_type_map[al_2]));

#define HYDRO 1


int cAtomicResidueMapper::initResidueAtoms() {
	if (nbType == 169) {
		return 2 + initResidueAtoms167();
		
	}
	if (nbType == 167) {
		return initResidueAtoms167();
	}
	else if (nbType == 4) {
		return initResidueAtoms4();
	}
	else if (nbType == 20) {
		return initResidueAtoms20();
	}
	else if (nbType == 21) {
		return initResidueAtoms21();
	}
	else {
		std::cout << "ERROR bad number of types : " << nbType << std::endl;
		exit(1);
	}
}

int cAtomicResidueMapper::initResidueAtoms4() {
	int n = 0;
	std::map<std::string, int> name_type_map;
	DEFINE_ATOMS("C", "N", "O", "S");
	//carbons
	DEFINE_ALIASE("CA" , "C");
	DEFINE_ALIASE("CB" , "C");
	DEFINE_ALIASE("CG" , "C");
	DEFINE_ALIASE("CG1", "C");
	DEFINE_ALIASE("CG2", "C");
	DEFINE_ALIASE("CD" , "C");
	DEFINE_ALIASE("CD1", "C");
	DEFINE_ALIASE("CD2", "C");
	DEFINE_ALIASE("CE" , "C");
	DEFINE_ALIASE("CE1", "C");
	DEFINE_ALIASE("CE2", "C");
	DEFINE_ALIASE("CH" , "C");
	DEFINE_ALIASE("CH1", "C");
	DEFINE_ALIASE("CH2", "C");
	DEFINE_ALIASE("CZ" , "C");
	DEFINE_ALIASE("CZ1", "C");
	DEFINE_ALIASE("CZ2", "C");
	DEFINE_ALIASE("CZ3", "C");


	//nitrogens	   N
	DEFINE_ALIASE("NA" , "N");
	DEFINE_ALIASE("NB" , "N");
	DEFINE_ALIASE("NG" , "N");
	DEFINE_ALIASE("NG1", "N");
	DEFINE_ALIASE("NG2", "N");
	DEFINE_ALIASE("ND" , "N");
	DEFINE_ALIASE("ND1", "N");
	DEFINE_ALIASE("ND2", "N");
	DEFINE_ALIASE("NE" , "N");
	DEFINE_ALIASE("NE1", "N");
	DEFINE_ALIASE("NE2", "N");
	DEFINE_ALIASE("NH" , "N");
	DEFINE_ALIASE("NH1", "N");
	DEFINE_ALIASE("NH2", "N");
	DEFINE_ALIASE("NZ" , "N");
	DEFINE_ALIASE("NZ1", "N");
	DEFINE_ALIASE("NZ2", "N");
	DEFINE_ALIASE("NZ3", "N");

	//oxygens	   O
	DEFINE_ALIASE("OA", "O");
	DEFINE_ALIASE("OB", "O");
	DEFINE_ALIASE("OG", "O");
	DEFINE_ALIASE("OG1", "O");
	DEFINE_ALIASE("OG2", "O");
	DEFINE_ALIASE("OD", "O");
	DEFINE_ALIASE("OD1", "O");
	DEFINE_ALIASE("OD2", "O");
	DEFINE_ALIASE("OE", "O");
	DEFINE_ALIASE("OE1", "O");
	DEFINE_ALIASE("OE2", "O");
	DEFINE_ALIASE("OH", "O");
	DEFINE_ALIASE("OH1", "O");
	DEFINE_ALIASE("OH2", "O");
	DEFINE_ALIASE("OZ", "O");
	DEFINE_ALIASE("OZ1", "O");
	DEFINE_ALIASE("OZ2", "O");
	DEFINE_ALIASE("OZ3", "O");
	DEFINE_ALIASE("OXT", "O");

	//sulfur	   O
	DEFINE_ALIASE("SG", "S");
	DEFINE_ALIASE("SD", "S");

	for (int c = 0; c < 23; c++) {
		atomTypes.insert(std::pair<cAminoResidue::eAminoType, std::map<std::string, int>>((cAminoResidue::eAminoType) c, name_type_map));
	}
	return n;
}

int cAtomicResidueMapper::initResidueAtoms21() {
	int n = 0;
	cAminoResidue::eAminoType cr;
	std::map<std::string, int> name_type_map;
	bool hydrogens = _hydrogens;

	cr = cAminoResidue::eAminoType::ALA;
	name_type_map.insert(std::pair<std::string, int>("CB", 0));

	atomTypes.insert(std::pair<cAminoResidue::eAminoType, std::map<std::string, int>>(cr, name_type_map));
	name_type_map.clear();

	cr = cAminoResidue::eAminoType::ARG;
	name_type_map.insert(std::pair<std::string, int>("CB", 1));

	atomTypes.insert(std::pair<cAminoResidue::eAminoType, std::map<std::string, int>>(cr, name_type_map));
	name_type_map.clear();

	cr = cAminoResidue::eAminoType::ASN;
	name_type_map.insert(std::pair<std::string, int>("CB", 2));

	atomTypes.insert(std::pair<cAminoResidue::eAminoType, std::map<std::string, int>>(cr, name_type_map));
	name_type_map.clear();

	cr = cAminoResidue::eAminoType::ASP;
	name_type_map.insert(std::pair<std::string, int>("CB", 3));

	atomTypes.insert(std::pair<cAminoResidue::eAminoType, std::map<std::string, int>>(cr, name_type_map));
	name_type_map.clear();

	cr = cAminoResidue::eAminoType::CYS;
	name_type_map.insert(std::pair<std::string, int>("CB", 4));

	atomTypes.insert(std::pair<cAminoResidue::eAminoType, std::map<std::string, int>>(cr, name_type_map));
	name_type_map.clear();

	cr = cAminoResidue::eAminoType::GLN;
	name_type_map.insert(std::pair<std::string, int>("CB", 5));

	atomTypes.insert(std::pair<cAminoResidue::eAminoType, std::map<std::string, int>>(cr, name_type_map));
	name_type_map.clear();

	cr = cAminoResidue::eAminoType::GLU;
	name_type_map.insert(std::pair<std::string, int>("CB", 6));

	atomTypes.insert(std::pair<cAminoResidue::eAminoType, std::map<std::string, int>>(cr, name_type_map));
	name_type_map.clear();

	cr = cAminoResidue::eAminoType::GLY;
	name_type_map.insert(std::pair<std::string, int>("CA", 7));

	atomTypes.insert(std::pair<cAminoResidue::eAminoType, std::map<std::string, int>>(cr, name_type_map));
	name_type_map.clear();

	cr = cAminoResidue::eAminoType::HIS;
	name_type_map.insert(std::pair<std::string, int>("CB", 8));

	atomTypes.insert(std::pair<cAminoResidue::eAminoType, std::map<std::string, int>>(cr, name_type_map));
	name_type_map.clear();

	cr = cAminoResidue::eAminoType::ILE;
	name_type_map.insert(std::pair<std::string, int>("CB", 9));

	atomTypes.insert(std::pair<cAminoResidue::eAminoType, std::map<std::string, int>>(cr, name_type_map));
	name_type_map.clear();

	cr = cAminoResidue::eAminoType::LEU;
	name_type_map.insert(std::pair<std::string, int>("CB", 10));

	atomTypes.insert(std::pair<cAminoResidue::eAminoType, std::map<std::string, int>>(cr, name_type_map));
	name_type_map.clear();

	cr = cAminoResidue::eAminoType::LYS;
	name_type_map.insert(std::pair<std::string, int>("CB", 11));

	atomTypes.insert(std::pair<cAminoResidue::eAminoType, std::map<std::string, int>>(cr, name_type_map));
	name_type_map.clear();

	cr = cAminoResidue::eAminoType::MET;
	name_type_map.insert(std::pair<std::string, int>("CB", 12));
	atomTypes.insert(std::pair<cAminoResidue::eAminoType, std::map<std::string, int>>(cr, name_type_map));
	cr = cAminoResidue::eAminoType::MSE;
	atomTypes.insert(std::pair<cAminoResidue::eAminoType, std::map<std::string, int>>(cr, name_type_map));
	name_type_map.clear();

	cr = cAminoResidue::eAminoType::PHE;
	name_type_map.insert(std::pair<std::string, int>("CB", 13));
	atomTypes.insert(std::pair<cAminoResidue::eAminoType, std::map<std::string, int>>(cr, name_type_map));
	name_type_map.clear();

	cr = cAminoResidue::eAminoType::PRO;
	name_type_map.insert(std::pair<std::string, int>("CB", 14));

	atomTypes.insert(std::pair<cAminoResidue::eAminoType, std::map<std::string, int>>(cr, name_type_map));
	name_type_map.clear();

	cr = cAminoResidue::eAminoType::SER;
	name_type_map.insert(std::pair<std::string, int>("CB", 15));

	atomTypes.insert(std::pair<cAminoResidue::eAminoType, std::map<std::string, int>>(cr, name_type_map));
	name_type_map.clear();

	cr = cAminoResidue::eAminoType::THR;
	name_type_map.insert(std::pair<std::string, int>("CB", 16));

	atomTypes.insert(std::pair<cAminoResidue::eAminoType, std::map<std::string, int>>(cr, name_type_map));
	name_type_map.clear();

	cr = cAminoResidue::eAminoType::TRP;
	name_type_map.insert(std::pair<std::string, int>("CB", 17));
	DEFINE_ALIASE("O1", "O");

	atomTypes.insert(std::pair<cAminoResidue::eAminoType, std::map<std::string, int>>(cr, name_type_map));
	name_type_map.clear();

	cr = cAminoResidue::eAminoType::TYR;
	name_type_map.insert(std::pair<std::string, int>("CB", 18));

	atomTypes.insert(std::pair<cAminoResidue::eAminoType, std::map<std::string, int>>(cr, name_type_map));
	name_type_map.clear();

	cr = cAminoResidue::eAminoType::VAL;
	name_type_map.insert(std::pair<std::string, int>("CB", 19));
	atomTypes.insert(std::pair<cAminoResidue::eAminoType, std::map<std::string, int>>(cr, name_type_map));

	name_type_map.clear();


	cr = cAminoResidue::eAminoType::AMINO_TYPE_END;
	name_type_map.insert(std::pair<std::string, int>("CB", 20));

	atomTypes.insert(std::pair<cAminoResidue::eAminoType, std::map<std::string, int>>(cr, name_type_map));
	return 21;
}


int cAtomicResidueMapper::initResidueAtoms20() {
	int n = 0;
	cAminoResidue::eAminoType cr;
	std::map<std::string, int> name_type_map;
	bool hydrogens = _hydrogens;

	cr = cAminoResidue::eAminoType::ALA;
	name_type_map.insert(std::pair<std::string, int>("C", 3));
	name_type_map.insert(std::pair<std::string, int>("N", 9));
	name_type_map.insert(std::pair<std::string, int>("O", 15));
	name_type_map.insert(std::pair<std::string, int>("CA", 7));
	name_type_map.insert(std::pair<std::string, int>("CB", 6));
	DEFINE_ALIASE("O1", "O");

	atomTypes.insert(std::pair<cAminoResidue::eAminoType, std::map<std::string, int>>(cr, name_type_map));
	name_type_map.clear();

	cr = cAminoResidue::eAminoType::ARG;
	name_type_map.insert(std::pair<std::string, int>("C", 3));
	name_type_map.insert(std::pair<std::string, int>("N", 9));
	name_type_map.insert(std::pair<std::string, int>("O", 15));
	name_type_map.insert(std::pair<std::string, int>("CA", 7));
	name_type_map.insert(std::pair<std::string, int>("CB", 6));
	name_type_map.insert(std::pair<std::string, int>("CG", 6));
	name_type_map.insert(std::pair<std::string, int>("CD", 8));
	name_type_map.insert(std::pair<std::string, int>("CZ", 1));
	name_type_map.insert(std::pair<std::string, int>("NE", 13));
	name_type_map.insert(std::pair<std::string, int>("NH1", 12));
	name_type_map.insert(std::pair<std::string, int>("NH2", 12));
	DEFINE_ALIASE("O1", "O");

	atomTypes.insert(std::pair<cAminoResidue::eAminoType, std::map<std::string, int>>(cr, name_type_map));
	name_type_map.clear();

	cr = cAminoResidue::eAminoType::ASN;
	name_type_map.insert(std::pair<std::string, int>("C", 3));
	name_type_map.insert(std::pair<std::string, int>("N", 9));
	name_type_map.insert(std::pair<std::string, int>("O", 15));
	name_type_map.insert(std::pair<std::string, int>("CA", 7));
	name_type_map.insert(std::pair<std::string, int>("CB", 6));
	name_type_map.insert(std::pair<std::string, int>("CG", 4));
	name_type_map.insert(std::pair<std::string, int>("OD1", 16));
	name_type_map.insert(std::pair<std::string, int>("ND2", 10));
	DEFINE_ALIASE("O1", "O");

	atomTypes.insert(std::pair<cAminoResidue::eAminoType, std::map<std::string, int>>(cr, name_type_map));
	name_type_map.clear();

	cr = cAminoResidue::eAminoType::ASP;
	name_type_map.insert(std::pair<std::string, int>("C", 3));
	name_type_map.insert(std::pair<std::string, int>("N", 9));
	name_type_map.insert(std::pair<std::string, int>("O", 15));
	name_type_map.insert(std::pair<std::string, int>("CA", 7));
	name_type_map.insert(std::pair<std::string, int>("CB", 6));
	name_type_map.insert(std::pair<std::string, int>("CG", 2));
	name_type_map.insert(std::pair<std::string, int>("OD1", 18));
	name_type_map.insert(std::pair<std::string, int>("OD2", 18));
	DEFINE_ALIASE("O1", "O");

	atomTypes.insert(std::pair<cAminoResidue::eAminoType, std::map<std::string, int>>(cr, name_type_map));
	name_type_map.clear();

	cr = cAminoResidue::eAminoType::CYS;
	name_type_map.insert(std::pair<std::string, int>("C", 3));
	name_type_map.insert(std::pair<std::string, int>("N", 9));
	name_type_map.insert(std::pair<std::string, int>("O", 15));
	name_type_map.insert(std::pair<std::string, int>("CA", 7));
	name_type_map.insert(std::pair<std::string, int>("CB", 8));
	name_type_map.insert(std::pair<std::string, int>("SG", 19));
	DEFINE_ALIASE("O1", "O");

	atomTypes.insert(std::pair<cAminoResidue::eAminoType, std::map<std::string, int>>(cr, name_type_map));
	name_type_map.clear();

	cr = cAminoResidue::eAminoType::GLN;
	name_type_map.insert(std::pair<std::string, int>("C", 3));
	name_type_map.insert(std::pair<std::string, int>("N", 9));
	name_type_map.insert(std::pair<std::string, int>("O", 15));
	name_type_map.insert(std::pair<std::string, int>("CA", 7));
	name_type_map.insert(std::pair<std::string, int>("CB", 6));
	name_type_map.insert(std::pair<std::string, int>("CG", 6));
	name_type_map.insert(std::pair<std::string, int>("CD", 4));
	name_type_map.insert(std::pair<std::string, int>("OE1", 16));
	name_type_map.insert(std::pair<std::string, int>("NE2", 10));
	DEFINE_ALIASE("O1", "O");

	atomTypes.insert(std::pair<cAminoResidue::eAminoType, std::map<std::string, int>>(cr, name_type_map));
	name_type_map.clear();

	cr = cAminoResidue::eAminoType::GLU;
	name_type_map.insert(std::pair<std::string, int>("C", 3));
	name_type_map.insert(std::pair<std::string, int>("N", 9));
	name_type_map.insert(std::pair<std::string, int>("O", 15));
	name_type_map.insert(std::pair<std::string, int>("CA", 7));
	name_type_map.insert(std::pair<std::string, int>("CB", 6));
	name_type_map.insert(std::pair<std::string, int>("CG", 6));
	name_type_map.insert(std::pair<std::string, int>("CD", 2));
	name_type_map.insert(std::pair<std::string, int>("OE1", 18));
	name_type_map.insert(std::pair<std::string, int>("OE2", 18));
	DEFINE_ALIASE("O1", "O");
	
	atomTypes.insert(std::pair<cAminoResidue::eAminoType, std::map<std::string, int>>(cr, name_type_map));
	name_type_map.clear();

	cr = cAminoResidue::eAminoType::GLY;
	name_type_map.insert(std::pair<std::string, int>("C", 3));
	name_type_map.insert(std::pair<std::string, int>("N", 9));
	name_type_map.insert(std::pair<std::string, int>("O", 15));
	name_type_map.insert(std::pair<std::string, int>("CA", 7));
	DEFINE_ALIASE("O1", "O");

	atomTypes.insert(std::pair<cAminoResidue::eAminoType, std::map<std::string, int>>(cr, name_type_map));
	name_type_map.clear();

	cr = cAminoResidue::eAminoType::HIS;
	name_type_map.insert(std::pair<std::string, int>("C", 3));
	name_type_map.insert(std::pair<std::string, int>("N", 9));
	name_type_map.insert(std::pair<std::string, int>("O", 15));
	name_type_map.insert(std::pair<std::string, int>("CA", 7));
	name_type_map.insert(std::pair<std::string, int>("CB", 6));
	name_type_map.insert(std::pair<std::string, int>("CG", 5));
	name_type_map.insert(std::pair<std::string, int>("ND1", 11));
	name_type_map.insert(std::pair<std::string, int>("CD2", 5));
	name_type_map.insert(std::pair<std::string, int>("CE1", 5));
	name_type_map.insert(std::pair<std::string, int>("NE2", 11));
	DEFINE_ALIASE("O1", "O");

	atomTypes.insert(std::pair<cAminoResidue::eAminoType, std::map<std::string, int>>(cr, name_type_map));
	name_type_map.clear();

	cr = cAminoResidue::eAminoType::ILE;
	name_type_map.insert(std::pair<std::string, int>("C", 3));
	name_type_map.insert(std::pair<std::string, int>("N", 9));
	name_type_map.insert(std::pair<std::string, int>("O", 15));
	name_type_map.insert(std::pair<std::string, int>("CA", 7));
	name_type_map.insert(std::pair<std::string, int>("CB", 6));
	name_type_map.insert(std::pair<std::string, int>("CG1", 6));
	name_type_map.insert(std::pair<std::string, int>("CG2", 6));
	name_type_map.insert(std::pair<std::string, int>("CD1", 6));
	DEFINE_ALIASE("O1", "O");
	DEFINE_ALIASE("CD", "CD1");
	
	atomTypes.insert(std::pair<cAminoResidue::eAminoType, std::map<std::string, int>>(cr, name_type_map));
	name_type_map.clear();

	cr = cAminoResidue::eAminoType::LEU;
	name_type_map.insert(std::pair<std::string, int>("C", 3));
	name_type_map.insert(std::pair<std::string, int>("N", 9));
	name_type_map.insert(std::pair<std::string, int>("O", 15));
	name_type_map.insert(std::pair<std::string, int>("CA", 7));
	name_type_map.insert(std::pair<std::string, int>("CB", 6));
	name_type_map.insert(std::pair<std::string, int>("CG", 6));
	name_type_map.insert(std::pair<std::string, int>("CD1", 6));
	name_type_map.insert(std::pair<std::string, int>("CD2", 6));
	DEFINE_ALIASE("O1", "O");

	atomTypes.insert(std::pair<cAminoResidue::eAminoType, std::map<std::string, int>>(cr, name_type_map));
	name_type_map.clear();

	cr = cAminoResidue::eAminoType::LYS;
	name_type_map.insert(std::pair<std::string, int>("C", 3));
	name_type_map.insert(std::pair<std::string, int>("N", 9));
	name_type_map.insert(std::pair<std::string, int>("O", 15));
	name_type_map.insert(std::pair<std::string, int>("CA", 7));
	name_type_map.insert(std::pair<std::string, int>("CB", 6));
	name_type_map.insert(std::pair<std::string, int>("CG", 6));
	name_type_map.insert(std::pair<std::string, int>("CD", 6));
	name_type_map.insert(std::pair<std::string, int>("CE", 8));
	name_type_map.insert(std::pair<std::string, int>("NZ", 14));
	DEFINE_ALIASE("O1", "O");

	atomTypes.insert(std::pair<cAminoResidue::eAminoType, std::map<std::string, int>>(cr, name_type_map));
	name_type_map.clear();

	cr = cAminoResidue::eAminoType::MET;
	name_type_map.insert(std::pair<std::string, int>("C", 3));
	name_type_map.insert(std::pair<std::string, int>("N", 9));
	name_type_map.insert(std::pair<std::string, int>("O", 15));
	name_type_map.insert(std::pair<std::string, int>("CA", 7));
	name_type_map.insert(std::pair<std::string, int>("CB", 6));
	name_type_map.insert(std::pair<std::string, int>("CG", 8));
	name_type_map.insert(std::pair<std::string, int>("SD", 0));
	name_type_map.insert(std::pair<std::string, int>("CE", 8));
	DEFINE_ALIASE("O1", "O");
	DEFINE_ALIASE("SE", "SD");
	atomTypes.insert(std::pair<cAminoResidue::eAminoType, std::map<std::string, int>>(cr, name_type_map));
	cr = cAminoResidue::eAminoType::MSE;
	atomTypes.insert(std::pair<cAminoResidue::eAminoType, std::map<std::string, int>>(cr, name_type_map));
	name_type_map.clear();

	cr = cAminoResidue::eAminoType::PHE;
	name_type_map.insert(std::pair<std::string, int>("C", 3));
	name_type_map.insert(std::pair<std::string, int>("N", 9));
	name_type_map.insert(std::pair<std::string, int>("O", 15));
	name_type_map.insert(std::pair<std::string, int>("CA", 7));
	name_type_map.insert(std::pair<std::string, int>("CB", 6));
	name_type_map.insert(std::pair<std::string, int>("CG", 5));
	name_type_map.insert(std::pair<std::string, int>("CD1", 5));
	name_type_map.insert(std::pair<std::string, int>("CD2", 5));
	name_type_map.insert(std::pair<std::string, int>("CE1", 5));
	name_type_map.insert(std::pair<std::string, int>("CE2", 5));
	name_type_map.insert(std::pair<std::string, int>("CZ", 5));
	DEFINE_ALIASE("O1", "O");
	atomTypes.insert(std::pair<cAminoResidue::eAminoType, std::map<std::string, int>>(cr, name_type_map));
	name_type_map.clear();

	cr = cAminoResidue::eAminoType::PRO;
	name_type_map.insert(std::pair<std::string, int>("C", 3));
	name_type_map.insert(std::pair<std::string, int>("N", 9));
	name_type_map.insert(std::pair<std::string, int>("O", 15));
	name_type_map.insert(std::pair<std::string, int>("CA", 7));
	name_type_map.insert(std::pair<std::string, int>("CB", 6));
	name_type_map.insert(std::pair<std::string, int>("CG", 6));
	name_type_map.insert(std::pair<std::string, int>("CD", 8));
	DEFINE_ALIASE("O1", "O");

	atomTypes.insert(std::pair<cAminoResidue::eAminoType, std::map<std::string, int>>(cr, name_type_map));
	name_type_map.clear();

	cr = cAminoResidue::eAminoType::SER;
	name_type_map.insert(std::pair<std::string, int>("C", 3));
	name_type_map.insert(std::pair<std::string, int>("N", 9));
	name_type_map.insert(std::pair<std::string, int>("O", 15));
	name_type_map.insert(std::pair<std::string, int>("CA", 7));
	name_type_map.insert(std::pair<std::string, int>("CB", 6));
	name_type_map.insert(std::pair<std::string, int>("OG", 17));
	DEFINE_ALIASE("O1", "O");

	atomTypes.insert(std::pair<cAminoResidue::eAminoType, std::map<std::string, int>>(cr, name_type_map));
	name_type_map.clear();

	cr = cAminoResidue::eAminoType::THR;
	name_type_map.insert(std::pair<std::string, int>("C", 3));
	name_type_map.insert(std::pair<std::string, int>("N", 9));
	name_type_map.insert(std::pair<std::string, int>("O", 15));
	name_type_map.insert(std::pair<std::string, int>("CA", 7));
	name_type_map.insert(std::pair<std::string, int>("CB", 8));
	name_type_map.insert(std::pair<std::string, int>("OG1", 17));
	name_type_map.insert(std::pair<std::string, int>("CG2", 6));
	DEFINE_ALIASE("O1", "O");
	atomTypes.insert(std::pair<cAminoResidue::eAminoType, std::map<std::string, int>>(cr, name_type_map));
	name_type_map.clear();

	cr = cAminoResidue::eAminoType::TRP;
	//DEFINE_ATOMS("C", "N", "O", "CA", "CB", "CD1", "CD2", "CE2", "CE3", "NE1", "CG", "CH2", "CZ2", "CZ3");
	name_type_map.insert(std::pair<std::string, int>("C", 3));
	name_type_map.insert(std::pair<std::string, int>("N", 9));
	name_type_map.insert(std::pair<std::string, int>("O", 15));
	name_type_map.insert(std::pair<std::string, int>("CA", 7));
	name_type_map.insert(std::pair<std::string, int>("CB", 6));
	name_type_map.insert(std::pair<std::string, int>("CG", 5));
	name_type_map.insert(std::pair<std::string, int>("CD1", 5));
	name_type_map.insert(std::pair<std::string, int>("CD2", 5));
	name_type_map.insert(std::pair<std::string, int>("CE2", 5));
	name_type_map.insert(std::pair<std::string, int>("CE3", 5));
	name_type_map.insert(std::pair<std::string, int>("NE1", 11));
	name_type_map.insert(std::pair<std::string, int>("CH2", 5));
	name_type_map.insert(std::pair<std::string, int>("CZ2", 5));
	name_type_map.insert(std::pair<std::string, int>("CZ3", 5));
	DEFINE_ALIASE("O1", "O");
	
	atomTypes.insert(std::pair<cAminoResidue::eAminoType, std::map<std::string, int>>(cr, name_type_map));
	name_type_map.clear();

	cr = cAminoResidue::eAminoType::TYR;
	//DEFINE_ATOMS("C", "N", "O", "CA", "CB", "CD1", "CD2", "CE1", "CE2", "CG", "OH", "CZ");
	name_type_map.insert(std::pair<std::string, int>("C", 3));
	name_type_map.insert(std::pair<std::string, int>("N", 9));
	name_type_map.insert(std::pair<std::string, int>("O", 15));
	name_type_map.insert(std::pair<std::string, int>("CA", 7));
	name_type_map.insert(std::pair<std::string, int>("CB", 6));
	name_type_map.insert(std::pair<std::string, int>("CG", 5));
	name_type_map.insert(std::pair<std::string, int>("CD1", 5));
	name_type_map.insert(std::pair<std::string, int>("CD2", 5));
	name_type_map.insert(std::pair<std::string, int>("CE1", 5));
	name_type_map.insert(std::pair<std::string, int>("CE2", 5));
	name_type_map.insert(std::pair<std::string, int>("OH", 17));
	name_type_map.insert(std::pair<std::string, int>("CZ", 5));
	DEFINE_ALIASE("O1", "O");

	atomTypes.insert(std::pair<cAminoResidue::eAminoType, std::map<std::string, int>>(cr, name_type_map));
	name_type_map.clear();

	cr = cAminoResidue::eAminoType::VAL;
	//DEFINE_ATOMS("C", "N", "O", "CA", "CB", "CG1", "CG2");
	name_type_map.insert(std::pair<std::string, int>("C", 3));
	name_type_map.insert(std::pair<std::string, int>("N", 9));
	name_type_map.insert(std::pair<std::string, int>("O", 15));
	name_type_map.insert(std::pair<std::string, int>("CA", 7));
	name_type_map.insert(std::pair<std::string, int>("CB", 6));
	name_type_map.insert(std::pair<std::string, int>("CG1", 6));
	name_type_map.insert(std::pair<std::string, int>("CG2", 6));
	DEFINE_ALIASE("O1", "O");
	atomTypes.insert(std::pair<cAminoResidue::eAminoType, std::map<std::string, int>>(cr, name_type_map));

	name_type_map.clear();


	cr = cAminoResidue::eAminoType::AMINO_TYPE_END;

	atomTypes.insert(std::pair<cAminoResidue::eAminoType, std::map<std::string, int>>(cr, name_type_map));
	return 20;
}

int cAtomicResidueMapper::initResidueAtoms167() {
	int n = 0;
	cAminoResidue::eAminoType cr;
	std::map<std::string, int> name_type_map;
	bool hydrogens = _hydrogens;

	cr = cAminoResidue::eAminoType::ALA;
	DEFINE_ATOMS("C", "N", "O", "CA", "CB");
	if (hydrogens) {
		DEFINE_ATOMS("H", "HA", "HB3", "HB2", "HB1");
		DEFINE_ALIASE("HA2", "HA");
	}
	DEFINE_ALIASE("O1", "O");

	atomTypes.insert(std::pair<cAminoResidue::eAminoType, std::map<std::string, int>>(cr, name_type_map));
	name_type_map.clear();

	cr = cAminoResidue::eAminoType::ARG;
	DEFINE_ATOMS("C", "N", "O", "CA", "CB", "CD", "NE", "CG", "NH1", "NH2", "CZ");
	DEFINE_ALIASE("O1", "O");
	if (hydrogens) {
		DEFINE_ATOMS("H", "HA", "HB2", "HB1", "HD2", "HD1", "HE", "HG2", "HG1", "HH11", "HH21", "HH12", "HH22");
		DEFINE_ALIASE("HA2", "HA");
		DEFINE_ALIASE("HB3", "HB1");
		DEFINE_ALIASE("HG3", "HG1");
		DEFINE_ALIASE("HD3", "HD1");
	}
	atomTypes.insert(std::pair<cAminoResidue::eAminoType, std::map<std::string, int>>(cr, name_type_map));
	name_type_map.clear();

	cr = cAminoResidue::eAminoType::ASN;
	DEFINE_ATOMS("C", "N", "O", "CA", "CB", "ND2", "OD1", "CG");
	DEFINE_ALIASE("O1", "O");
	if (hydrogens) {
		DEFINE_ATOMS("H", "HA", "HB2", "HB1", "HD22", "HD21");
		DEFINE_ALIASE("HA2", "HA");
		DEFINE_ALIASE("HB3", "HB1");
	}
	atomTypes.insert(std::pair<cAminoResidue::eAminoType, std::map<std::string, int>>(cr, name_type_map));
	name_type_map.clear();

	cr = cAminoResidue::eAminoType::ASP;
	DEFINE_ATOMS("C", "N", "O", "CA", "CB", "OD1", "OD2", "CG");
	if (hydrogens) {
		DEFINE_ATOMS("H", "HA", "HB2", "HB1");
		DEFINE_ALIASE("HA2", "HA");
		DEFINE_ALIASE("HB3", "HB1");
	}
	DEFINE_ALIASE("O1", "O");
	atomTypes.insert(std::pair<cAminoResidue::eAminoType, std::map<std::string, int>>(cr, name_type_map));
	name_type_map.clear();

	cr = cAminoResidue::eAminoType::CYS;
	DEFINE_ATOMS("C", "N", "O", "CA", "CB", "SG");
	DEFINE_ALIASE("O1", "O");
	if (hydrogens) {
		DEFINE_ATOMS("H", "HA", "HB2", "HB1", "HG");
		DEFINE_ALIASE("HA2", "HA");
		DEFINE_ALIASE("HG1", "HG");
		DEFINE_ALIASE("HB3", "HB1");
	}
	atomTypes.insert(std::pair<cAminoResidue::eAminoType, std::map<std::string, int>>(cr, name_type_map));
	name_type_map.clear();

	cr = cAminoResidue::eAminoType::GLN;
	DEFINE_ATOMS("C", "N", "O", "CA", "CB", "CD", "NE2", "OE1", "CG");
	DEFINE_ALIASE("O1", "O");
	if (hydrogens) {
		DEFINE_ATOMS("H", "HA", "HB2", "HB1", "HG2", "HG1", "HE22", "HE21");
		DEFINE_ALIASE("HA2", "HA");
		DEFINE_ALIASE("HB3", "HB1");
		DEFINE_ALIASE("HG3", "HG1");
	}
	atomTypes.insert(std::pair<cAminoResidue::eAminoType, std::map<std::string, int>>(cr, name_type_map));
	name_type_map.clear();

	cr = cAminoResidue::eAminoType::GLU;
	DEFINE_ATOMS("C", "N", "O", "CA", "CB", "CD", "OE1", "OE2", "CG");
	DEFINE_ALIASE("O1", "O");
	if (hydrogens) {
		DEFINE_ATOMS("H", "HA", "HB2", "HB1", "HG2", "HG1");
		DEFINE_ALIASE("HA2", "HA");
		DEFINE_ALIASE("HB3", "HB1");
		DEFINE_ALIASE("HG3", "HG1");
	}
	atomTypes.insert(std::pair<cAminoResidue::eAminoType, std::map<std::string, int>>(cr, name_type_map));
	name_type_map.clear();

	cr = cAminoResidue::eAminoType::GLY;
	DEFINE_ATOMS("C", "N", "O", "CA");
	DEFINE_ALIASE("O1", "O");
	if (hydrogens) {
		DEFINE_ATOMS("H", "HA1", "HA2");
		DEFINE_ALIASE("HA3", "HA1");
		DEFINE_ALIASE("HA", "HA1");
	}
	atomTypes.insert(std::pair<cAminoResidue::eAminoType, std::map<std::string, int>>(cr, name_type_map));
	name_type_map.clear();

	cr = cAminoResidue::eAminoType::HIS;
	DEFINE_ATOMS("C", "N", "O", "CA", "CB", "CD2", "ND1", "CE1", "NE2", "CG");
	DEFINE_ALIASE("O1", "O");
	if (hydrogens) {
		DEFINE_ATOMS("H", "HA", "HB2", "HB1", "HD1", "HD2", "HE1", "HE2");
		DEFINE_ALIASE("HA2", "HA");
		DEFINE_ALIASE("HB3", "HB1");
		DEFINE_ALIASE("HD3", "HD1");
		DEFINE_ALIASE("HE3", "HE1");
	}
	atomTypes.insert(std::pair<cAminoResidue::eAminoType, std::map<std::string, int>>(cr, name_type_map));
	name_type_map.clear();

	cr = cAminoResidue::eAminoType::ILE;
	DEFINE_ATOMS("C", "N", "O", "CA", "CB", "CD1", "CG1", "CG2");
	DEFINE_ALIASE("O1", "O");
	DEFINE_ALIASE("CD", "CD1");
	if (hydrogens) {
		DEFINE_ATOMS("H", "HA", "HB", "HD13", "HD12", "HD11", "HG12", "HG23", "HG11", "HG22", "HG21");
		DEFINE_ALIASE("HA2", "HA");
		DEFINE_ALIASE("HG13", "HG11");
		DEFINE_ALIASE("HD1", "HD11");
		DEFINE_ALIASE("HD2", "HD12");
		DEFINE_ALIASE("HD3", "HD13");
	}
	atomTypes.insert(std::pair<cAminoResidue::eAminoType, std::map<std::string, int>>(cr, name_type_map));
	name_type_map.clear();

	cr = cAminoResidue::eAminoType::LEU;
	DEFINE_ATOMS("C", "N", "O", "CA", "CB", "CD1", "CD2", "CG");
	DEFINE_ALIASE("O1", "O");
	if (hydrogens) {
		DEFINE_ATOMS("H", "HA", "HB2", "HB1", "HD13", "HD23", "HD12", "HD22", "HD11", "HD21", "HG");
		DEFINE_ALIASE("HA2", "HA");
		DEFINE_ALIASE("HB3", "HB1");
	}
	atomTypes.insert(std::pair<cAminoResidue::eAminoType, std::map<std::string, int>>(cr, name_type_map));
	name_type_map.clear();

	cr = cAminoResidue::eAminoType::LYS;
	DEFINE_ATOMS("C", "N", "O", "CA", "CB", "CD", "CE", "CG", "NZ");
	DEFINE_ALIASE("O1", "O");
	if (hydrogens) {
		DEFINE_ATOMS("H", "HA", "HB2", "HB1", "HD2", "HD1", "HE2", "HE1", "HG2", "HG1", "HZ3", "HZ2", "HZ1");
		DEFINE_ALIASE("HA2", "HA");
		DEFINE_ALIASE("HB3", "HB1");
		DEFINE_ALIASE("HG3", "HG1");
		DEFINE_ALIASE("HD3", "HD1");
		DEFINE_ALIASE("HE3", "HE1");
	}
	atomTypes.insert(std::pair<cAminoResidue::eAminoType, std::map<std::string, int>>(cr, name_type_map));
	name_type_map.clear();

	cr = cAminoResidue::eAminoType::MET;
	DEFINE_ATOMS("C", "N", "O", "CA", "CB", "CE", "CG", "SD");
	DEFINE_ALIASE("O1", "O");
	DEFINE_ALIASE("SE", "SD");
	if (hydrogens) {
		DEFINE_ATOMS("H", "HA", "HB2", "HB1", "HE3", "HE2", "HE1", "HG2", "HG1");
		DEFINE_ALIASE("HA2", "HA");
		DEFINE_ALIASE("HG3", "HG1");
		DEFINE_ALIASE("HB3", "HB1");
	}
	atomTypes.insert(std::pair<cAminoResidue::eAminoType, std::map<std::string, int>>(cr, name_type_map));
	cr = cAminoResidue::eAminoType::MSE;
	atomTypes.insert(std::pair<cAminoResidue::eAminoType, std::map<std::string, int>>(cr, name_type_map));
	name_type_map.clear();

	cr = cAminoResidue::eAminoType::PHE;
	DEFINE_ATOMS("C", "N", "O", "CA", "CB", "CD1", "CD2", "CE1", "CE2", "CG", "CZ");
	DEFINE_ALIASE("O1", "O");
	if (hydrogens) {
		DEFINE_ATOMS("H", "HA", "HB2", "HB1", "HD1", "HD2", "HE1", "HE2", "HZ");
		DEFINE_ALIASE("HA2", "HA");
		DEFINE_ALIASE("HB3", "HB1");
	}
	atomTypes.insert(std::pair<cAminoResidue::eAminoType, std::map<std::string, int>>(cr, name_type_map));
	name_type_map.clear();

	cr = cAminoResidue::eAminoType::PRO;
	DEFINE_ATOMS("C", "N", "O", "CA", "CB", "CD", "CG");
	DEFINE_ALIASE("O1", "O");
	if (hydrogens) {
		DEFINE_ATOMS("HA", "HB2", "HB1", "HG2", "HG1", "HD2", "HD1");
		DEFINE_ALIASE("HA2", "HA");
		DEFINE_ALIASE("HB3", "HB1");
		DEFINE_ALIASE("HD3", "HD1");
		DEFINE_ALIASE("HG3", "HG1");
	}
	atomTypes.insert(std::pair<cAminoResidue::eAminoType, std::map<std::string, int>>(cr, name_type_map));
	name_type_map.clear();

	cr = cAminoResidue::eAminoType::SER;
	DEFINE_ATOMS("C", "N", "O", "CA", "CB", "OG");
	DEFINE_ALIASE("O1", "O");
	if (hydrogens) {
		DEFINE_ATOMS("H", "HA", "HB2", "HB1", "HG");
		DEFINE_ALIASE("HA2", "HA");
		DEFINE_ALIASE("HG1", "HG");
		DEFINE_ALIASE("HB3", "HB1");
	}
	atomTypes.insert(std::pair<cAminoResidue::eAminoType, std::map<std::string, int>>(cr, name_type_map));
	name_type_map.clear();

	cr = cAminoResidue::eAminoType::THR;
	DEFINE_ATOMS("C", "N", "O", "CA", "CB", "CG2", "OG1");
	DEFINE_ALIASE("O1", "O");
	if (hydrogens) {
		DEFINE_ATOMS("H", "HA", "HB", "HG1", "HG23", "HG22", "HG21");
		DEFINE_ALIASE("HA2", "HA");
	}
	atomTypes.insert(std::pair<cAminoResidue::eAminoType, std::map<std::string, int>>(cr, name_type_map));
	name_type_map.clear();

	cr = cAminoResidue::eAminoType::TRP;
	DEFINE_ATOMS("C", "N", "O", "CA", "CB", "CD1", "CD2", "CE2", "CE3", "NE1", "CG", "CH2", "CZ2", "CZ3");
	DEFINE_ALIASE("O1", "O");
	if (hydrogens) {
		DEFINE_ATOMS("H", "HA", "HB2", "HB1", "HD1", "HE1", "HE3", "HH2", "HZ2", "HZ3");
		DEFINE_ALIASE("HA2", "HA");
		DEFINE_ALIASE("HB3", "HB1");
		DEFINE_ALIASE("HZ1", "HZ3");
	}
	atomTypes.insert(std::pair<cAminoResidue::eAminoType, std::map<std::string, int>>(cr, name_type_map));
	name_type_map.clear();

	cr = cAminoResidue::eAminoType::TYR;
	DEFINE_ATOMS("C", "N", "O", "CA", "CB", "CD1", "CD2", "CE1", "CE2", "CG", "OH", "CZ");
	DEFINE_ALIASE("O1", "O");
	if (hydrogens) {
		DEFINE_ATOMS("H", "HA", "HB2", "HB1", "HD1", "HD2", "HE1", "HE2", "HH2");
		DEFINE_ALIASE("HA2", "HA");
		DEFINE_ALIASE("HB3", "HB1");
		DEFINE_ALIASE("HD3", "HD1");
		DEFINE_ALIASE("HE3", "HE1");
	}
	atomTypes.insert(std::pair<cAminoResidue::eAminoType, std::map<std::string, int>>(cr, name_type_map));
	name_type_map.clear();
	
	cr = cAminoResidue::eAminoType::VAL;
	DEFINE_ATOMS("C", "N", "O", "CA", "CB", "CG1", "CG2");
	DEFINE_ALIASE("O1", "O");
	if (hydrogens) {
		DEFINE_ATOMS("H", "HA", "HB", "HG13", "HG23", "HG12", "HG22", "HG11", "HG21");
		DEFINE_ALIASE("HA2", "HA");
	}
	atomTypes.insert(std::pair<cAminoResidue::eAminoType, std::map<std::string, int>>(cr, name_type_map));
	name_type_map.clear();


	cr = cAminoResidue::eAminoType::AMINO_TYPE_END;
	//DEFINE_ATOMS("OXT");
	//DEFINE_ALIASE("O1", "O");
	if (hydrogens) {
		DEFINE_ATOMS("H", "HA", "HB", "HG13", "HG23", "HG12", "HG22", "HG11", "HG21");
		DEFINE_ALIASE("HA2", "HA");
	}
	atomTypes.insert(std::pair<cAminoResidue::eAminoType, std::map<std::string, int>>(cr, name_type_map));

//std::cout << "Atom types initialized, " << n << " atom types created" << std::endl;

return n;
}


void cCoarseResidueMapper::setProtein(cProtein *protein) {
	cResidueMapper::setProtein(protein);

	currentResIt = protein->residues().begin();
	neighbourRes = NeighbourResidues(protein, sqrt(3)*0.5* gridSize *voxWidth + 3 * resSpread, 1);
}

void cResidueMapper::setProtein(cProtein *protein) {
	cProteinMapper::setProtein(protein);

	currentResIt = protein->residues().begin();


}

void cAtomicResidueMapper::setProtein(cProtein *protein) {
	cResidueMapper::setProtein(protein);


	std::vector<cAtom*> atoms;
	for (cAtom at : protein->atoms()) {
		atoms.push_back(new cAtom(at));
	}
	//TODO: fix memory leak here
	cAtom** PDBatoms = new cAtom*[atoms.size()];
	for (int i = 0; i < atoms.size(); i++) {
		PDBatoms[i] = atoms[i];
	}
	neighbourGrid.initPointer(PDBatoms, atoms.size(), NULL, sqrt(3)*0.5* gridSize *voxWidth + sqrt(log(255)) * resSpread);


}

cAtomicResidueMapper::~cAtomicResidueMapper()
{
	delete[] _map;
	delete[] _meta;
}

int cResidueMapper::getNbMaps(bool native, std::string scoreFilename)
{
	//Number of maps is the number of selected residues
	///////////////////////////////////////////
	return _resList.size();
	///////////////////////////////////////////

	int nMap = 0;
		int nRes = _protein->numResidues();

		auto residue = _protein->residues().begin();
		int resCode = -1;
		bool eof = false;
		std::string line;
		std::ifstream scoreFile;
		if (!scoreFilename.empty()) {
			scoreFile.open(scoreFilename, std::ios_base::in);
		}
		float score;
		for (int indexMap = 0; indexMap < nRes; indexMap++) {
			if (!scoreFilename.empty()) {
				while (resCode < (int)residue->seqNumber && !eof && scoreFile) {
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
			else if (resCode == residue->seqNumber) {
				scoreInt = 1000000 * score;
			}
			else {
				scoreInt = -1000000;
			}
			if (scoreInt >= 0) {
				nMap++;
			}
			++residue;
		}
	
	return nMap;
}
#define NB_THREADS 1

void cAtomicResidueMapper::runCurrentThread( int threadId){
	
	
	cVector3 middle = origin + (vX+vY+vZ) *(gridSize*voxWidth*0.5);
	std::set<cAtom * > neighbours;
	neighbourGrid.constructNeighbours(middle, neighbours);
#pragma omp parallel for
	for (auto atom : neighbours) {
		addAtom(atom);
	}
	
}

void cAtomicResidueMapper::runCurrent() {
	
	memset(_map, 0, gridSize*gridSize *gridSize *(nbType+6)); //+4 for SCR, +2 for R/L

	//std::cout<<"chainId\t"<<currentResIt->chainId<<std::endl;
	//std::cout<<"seqNumber\t"<<currentResIt->seqNumber<<std::endl;
	//std::cout<<"resName\t"<<currentResIt->resName<<std::endl;
	srand ( time(NULL) );
	int chosenNumber = (rand() % currentResIt->numAtoms());
	//std::cout<<'1'<<'\t'<<centerSphere[0]<<"\t"<<centerSphere[1]<<"\t"<<centerSphere[2]<<'\t'<<chosenNumber<<'\t'<<currentResIt->numAtoms()<<std::endl;		
	centerSphere = currentResIt->getAtom(chosenNumber)->getPosition();
	//std::cout<<'2'<<'\t'<<centerSphere[0]<<"\t"<<centerSphere[1]<<"\t"<<centerSphere[2]<<'\t'<<chosenNumber<<'\t'<<currentResIt->numAtoms()<<std::endl;		

	std::map<_pair,char>::iterator it = _resList.find({currentResIt->seqNumber,currentResIt->chainId});
	if(it != _resList.end())
	{
		_skipRes = false;
	}
	else
	{
		_skipRes = true;
		return;
	}

	setLocalFrame(&(*currentResIt), NULL, _orient);
	if (currentResIt->type() == cResidue::eType::AMINO) {

		std::vector<std::thread> threads;
		for (int j = 1; j < NB_THREADS; j++) {
			threads.emplace_back(&cAtomicResidueMapper::runCurrentThread, this, j);
		}
		runCurrentThread(0);
		cAminoResidue * aminoRes = (cAminoResidue*) &(*currentResIt);
		((int32_t*)_meta)[0] = aminoRes->seqNumber;
		_meta[1] = aminoRes->getResCode();
		_meta[2] = -1000000;

		for (int j = 0; j < NB_THREADS-1; j++) {
			threads[j].join();
		}
	}
}


void cAtomicResidueMapper::runGradient(float* grads, float* retypeMat, std::unordered_map<size_t, cVector3> &gradByAtom) {
	memset(_map, 0, gridSize*gridSize *gridSize *(nbType+6)); //+4 for SCR, +2 for R/L
	cAtom frameAtoms[5];
	setLocalFrame(&(*currentResIt), frameAtoms);
	
	cVector3 sumGrads(0,0,0);
	cVector3 rotGrad(0, 0, 0);
	if (currentResIt->type() == cResidue::eType::AMINO) {
		cAminoResidue * aminoRes = (cAminoResidue*) &(*currentResIt);
		cVector3 middle = origin + (vX + vY + vZ) *(gridSize*voxWidth*0.5);
		std::set<cAtom * > neighbours;
		neighbourGrid.constructNeighbours(middle, neighbours);
		for (auto atom : neighbours) {
			int atomType = getType(atom);
			cVector3 tmpGrad = addGradAtom(atom, grads, &(retypeMat[numRetype*atomType]));
			sumGrads += tmpGrad;
			rotGrad += atom->getPosition() ^ tmpGrad;
			gradByAtom[atom->atomSerial] += tmpGrad;
		}
		
	}

	cVector3 tmpGrad = -sumGrads;
	sumGrads += tmpGrad;
	rotGrad += frameAtoms[0].getPosition() ^ tmpGrad;
	gradByAtom[frameAtoms[0].atomSerial] += tmpGrad;

	//	rotation around local X axis
	double rotX = (worldToLocal * rotGrad)[0];
	double dy = (worldToLocal * (frameAtoms[4].getPosition() - frameAtoms[3].getPosition()))[1];
	double gradZdiff = rotX / dy;

	tmpGrad = vZ*gradZdiff;
	sumGrads += tmpGrad;
	rotGrad += frameAtoms[3].getPosition() ^ tmpGrad;
	gradByAtom[frameAtoms[3].atomSerial] += tmpGrad;

	tmpGrad = -vZ*gradZdiff;
	sumGrads += tmpGrad;
	rotGrad += frameAtoms[4].getPosition() ^ tmpGrad;
	gradByAtom[frameAtoms[4].atomSerial] += tmpGrad;


	//	rotation around local Y axis
	double rotY = (worldToLocal * rotGrad)[1];
	double dx = (worldToLocal * (frameAtoms[2].getPosition() - frameAtoms[1].getPosition()))[0];
	gradZdiff = rotY / dx;

	tmpGrad = -vZ*gradZdiff;
	sumGrads += tmpGrad;
	rotGrad += frameAtoms[1].getPosition() ^ tmpGrad;
	gradByAtom[frameAtoms[1].atomSerial] += tmpGrad;

	tmpGrad = vZ*gradZdiff;
	sumGrads += tmpGrad;
	rotGrad += frameAtoms[2].getPosition() ^ tmpGrad;
	gradByAtom[frameAtoms[2].atomSerial] += tmpGrad;


	//	rotation around local Z axis
	double rotZ = (worldToLocal * rotGrad)[2];
	double gradYdiff = rotZ / dx;

	tmpGrad = vY*gradYdiff;
	sumGrads += tmpGrad;
	rotGrad += frameAtoms[1].getPosition() ^ tmpGrad;
	gradByAtom[frameAtoms[1].atomSerial] += tmpGrad;

	tmpGrad = -vY*gradYdiff;
	sumGrads += tmpGrad;
	rotGrad += frameAtoms[2].getPosition() ^ tmpGrad;
	gradByAtom[frameAtoms[2].atomSerial] += tmpGrad;

	cVector3 sumGradAll(0, 0, 0);
	for (auto it = gradByAtom.begin(); it != gradByAtom.end(); ++it) {
		sumGradAll += it->second;
	}

	assert(sumGradAll.norm2() < 1e-15);
	return;

	
}


/*
void cAtomicResidueMapper::runGradient2(float* grads, float* retypeMat, std::unordered_map<size_t, cVector3> &gradByAtom) {
	memset(_map, 0, gridSize*gridSize *gridSize *nbType);
	setLocalFrame(&(*currentResIt));
	if (currentResIt->type() == cResidue::eType::AMINO) {
		cAminoResidue * aminoRes = (cAminoResidue*) &(*currentResIt);
		cAtom * nitro = aminoRes->getAtom("N");

		cAtom * carbon = aminoRes->getAtom("C");
		cAtom * carbonA = aminoRes->getAtom("CA");

		cVector3 gradCA = addGradAtom(carbonA, grads, &(retypeMat[numRetype*getType(carbonA)]));
		cVector3 gradC = addGradAtom(carbon, grads, &(retypeMat[numRetype*getType(carbon)]));
		cVector3 originGrad(0, 0, 0);
		cVector3 momGrad(0, 0, 0);
		cVector3 tmpGrad;
		for (cAtom atomIt : aminoRes->atoms()) {
			if (gradByAtom.find(atomIt.atomSerial) != gradByAtom.end()
				&& atomIt.atomSerial != carbon->atomSerial
				&& atomIt.atomSerial != carbonA->atomSerial
				&& atomIt.atomSerial != nitro->atomSerial) {
				int atomType = getType(&atomIt);
				tmpGrad = addGradAtom(&atomIt, grads, &(retypeMat[numRetype*atomType]));
				gradByAtom[atomIt.atomSerial] += tmpGrad;
				originGrad -= tmpGrad;
				momGrad += (atomIt.getPosition() - origin) ^ tmpGrad;
			}
		};
		if (neighbourRes.find(aminoRes) != neighbourRes.end()) {
			std::vector<cResidue*> neighbour_residues = neighbourRes.at(aminoRes);
			for (const auto &residue_neighbours : neighbour_residues) {
				if (&(*residue_neighbours) == &(*aminoRes)) {

				}
				for (cAtom atomIt : residue_neighbours->atoms()) {//->atoms().begin(); atomIt != residue_neighbours->atoms().end(); ++atomIt) {

					if (gradByAtom.find(atomIt.atomSerial) != gradByAtom.end()) {
						int atomType = getType(&atomIt);
						tmpGrad = addGradAtom(&atomIt, grads, &(retypeMat[numRetype*atomType]));
						gradByAtom[atomIt.atomSerial] += tmpGrad;
						originGrad -= tmpGrad;
						momGrad += (atomIt.getPosition() - origin) ^ tmpGrad;
					}
				};
			}
		}
	}
}
*/
bool cResidueMapper::setLocalFrame(cResidue* res, cAtom * frameAtoms, bool orient) {
	if (!res) return false;
	if (!res->getAtom("N")) return false;
	if (!res->getAtom("CA")) return false;
	if (!res->getAtom("C")) return false;
	origin = res->getAtom("N")->getPosition();
	if (!orient) {
		vX = cVector3(1, 0, 0);
		vY = cVector3(0, 1, 0);
		vZ = cVector3(0, 0, 1);
		//put the origin at the bottom corner of the map;
		origin -= voxWidth * gridSize * 0.5 * (vX + vY + vZ);
		worldToLocal = cMatrix33(vX, vY, vZ).getTranspose();
		localToWorld = cMatrix33(vX, vY, vZ);

		return true;
	}
	if (frameAtoms != NULL) {
		frameAtoms[0] = res->getAtom("N");
	}
	cVector3 minusX;
	if (res->getPrevious() != NULL && res->getPrevious()->getAtom("C") != NULL) {
		minusX = origin - res->getPrevious()->getAtom("C")->getPosition();
		if (frameAtoms != NULL) {
			frameAtoms[1] = res->getAtom("N");
			frameAtoms[2] = res->getPrevious()->getAtom("C");
		}
	}
	else {
		minusX = res->getAtom("C")->getPosition() - res->getAtom("CA")->getPosition();
		if (frameAtoms != NULL) {
			frameAtoms[1] = res->getAtom("C");
			frameAtoms[2] = res->getAtom("CA");
		}
	}
	vX = -minusX.normalizedVersion();
	cVector3 tmpY = res->getAtom("CA")->getPosition() - origin;
	if (frameAtoms != NULL) {
		frameAtoms[3] = res->getAtom("N");
		frameAtoms[4] = res->getAtom("CA");
	}
	tmpY -= vX * (tmpY | vX);
	vY = tmpY.normalizedVersion();
	vZ = vX^vY;
	// center the origin in the estimated amino acid center
	origin += 3.5 * vX;
	origin += 3.0 * vY;

	//put the origin at the bottom corner of the map;
	origin -= voxWidth * gridSize * 0.5 * (vX + vY + vZ);
	worldToLocal = cMatrix33(vX, vY, vZ).getTranspose();
	localToWorld = cMatrix33(vX, vY, vZ);
	//cVector3 p = worldToLocal * vX;
	//p = worldToLocal * vY;
	//p = worldToLocal * vZ;
	return true;
}

//this is the version for oriented graphs
bool cResidueMapper::setLocalFrameCa(cResidue* res, cAtom * frameAtoms, bool orient) {
    if (!res) return false;
    if (!res->getAtom("N")) return false;
    if (!res->getAtom("CA")) return false;
    if (!res->getAtom("C")) return false;
    origin = res->getAtom("CA")->getPosition();
    if (!orient) {
        vX = cVector3(1, 0, 0);
        vY = cVector3(0, 1, 0);
        vZ = cVector3(0, 0, 1);

        worldToLocal = cMatrix33(vX, vY, vZ).getTranspose();
        localToWorld = cMatrix33(vX, vY, vZ);

        return true;
    }
    if (frameAtoms != NULL) {
        frameAtoms[0] = res->getAtom("CA");
    }

    // vector X is from CA to N
    vX = res->getAtom("N")->getPosition() - origin;
    if (frameAtoms != NULL) {
        frameAtoms[1] = res->getAtom("CA");
        frameAtoms[2] = res->getAtom("N");
    }
    vX.normalize();

    // vector Y is from CA to C
    vY = res->getAtom("C")->getPosition() - origin;
    if (frameAtoms != NULL) {
        frameAtoms[3] = res->getAtom("CA");
        frameAtoms[4] = res->getAtom("C");
    }
    vY -= vX * (vY | vX);
    vY.normalize();

    //vZ is orthogonal to the XY plane
    vZ = vX^vY;

    worldToLocal = cMatrix33(vX, vY, vZ).getTranspose();
    localToWorld = cMatrix33(vX, vY, vZ);
    return true;
}

inline
double exp256(double x) {
	x = 1.0 + x / 256.0;
	x *= x; x *= x; x *= x; x *= x;
	x *= x; x *= x; x *= x; x *= x;
	return x;
}

int cAtomicResidueMapper::addAtom(cAtom* atom) {

	// Mask Carbon atoms from the current residue
	////////////////////////////////////////////////////////////////////////////
	//if (atom->name[0] == 'C' && atom->resSeqNumber == currentResIt->seqNumber) {
	//	std::cout<<"maskC\t"<<currentResIt->seqNumber<<"\t"<<atom->resSeqNumber<<"\t"<<atom->name[0]<<std::endl;
	//	return 0;
	//}
	// if (atom->name[0] == 'C') {
	// 	std::cout<<"nomaskC\t"<<currentResIt->seqNumber<<"\t"<<atom->resSeqNumber<<"\t"<<atom->name[0]<<std::endl;
	// }

	std::cout<<"mask\t"<<currentResIt->seqNumber<<"\t"<<atom->resSeqNumber<<"\t"<<atom->name[0]<<"\t"<<centerSphere[0]<<"\t"<<centerSphere[1]<<"\t"<<centerSphere[2]<<std::endl;		
	cVector3 positionCurrentAtom = atom->getPosition();
	cVector3 v_dist = positionCurrentAtom - centerSphere;
	float dist = sqrt(pow(v_dist[0],2) + pow(v_dist[1],2) + pow(v_dist[2],2));
	std::cout<<"DIST\t"<<dist<<std::endl;
	if (dist < 3)
	{
		return 0;
	}

	////////////////////////////////////////////////////////////////////////////

	if (atom->elementSymbol[0] != 'C' && atom->elementSymbol[0] != 'O'&&  atom->elementSymbol[0] != 'N' && atom->elementSymbol[0] != 'S') {
		if (atom->elementSymbol[0] != 'H') {
			//std::cout << "Ignored atom " << atom << std::endl;
		}
		return 0;
	}
	if (atom->residue->type() != cResidue::eType::AMINO) {
		return 0;
	}
	if ( abs((int)(atom->resSeqNumber - currentResIt->seqNumber)) <= _skipNeighb && abs((int)(atom->resSeqNumber - currentResIt->seqNumber))>0) {

		return 0;
	}
	int type = getType(atom);
	double atomSpread2 = atomSpread * atomSpread;
	if (type >= 0) {
		cVector3 localPos = (atom->getPosition() - origin);
		double lx = (localPos | vX), ly = (localPos | vY), lz = (localPos | vZ);
		float gridLength = gridSize * voxWidth;
		float cutoffDist = sqrt(log(255))*atomSpread;

		double cutoff = sqrt(log(255));
		if (lx > -cutoffDist && lx < gridLength + cutoffDist &&
			ly > -cutoffDist && ly < gridLength + cutoffDist &&
			lz > -cutoffDist && lz < gridLength + cutoffDist) {
			int ix = lx / voxWidth;
			int iy = ly / voxWidth;
			int iz = lz / voxWidth;
			for (int ixx = std::max(0., ix - cutoff); ixx < ix + cutoff+1 && ixx < gridSize; ixx++) {
				for (int iyy = std::max(0., iy - cutoff); iyy < iy + cutoff+1 && iyy < gridSize; iyy++) {
					for (int izz = std::max(0., iz - cutoff); izz < iz + cutoff+1 && izz < gridSize; izz++) {
						_map[((type*gridSize + izz)*gridSize + iyy )* gridSize + ixx] += 255 * exp256(-(
							(ixx * voxWidth - lx) * (ixx * voxWidth - lx) +
							(iyy * voxWidth - ly) * (iyy * voxWidth - ly) +
							(izz * voxWidth - lz) * (izz * voxWidth - lz)
							) / atomSpread2);
						
					}
				}
			}


			// SCR
			/////////////////////////////////////////////
			int scrType = 170;
			// std::cout<<"atominfo: "<<atom->resSeqNumber<<"\t"<<atom->chainId<<std::endl;
			std::map<_pair,char>::iterator it = _mapSCR.find({atom->resSeqNumber,atom->chainId});
			if(it != _mapSCR.end())
			{
				switch (it->second)
				{
				case 'S':
					scrType = 167;
					break;
				case 'C':
					scrType = 168;
					break;
				case 'R':
					scrType = 169;
					break;
				default:
					break;
				}
			}
			else{
				// std::cout<<"Notfound"<<std::endl;	
			}
			// std::cout<<"scrType: "<<scrType<<std::endl;
			
			for (int ixx = std::max(0., ix - cutoff); ixx < ix + cutoff+1 && ixx < gridSize; ixx++) {
				for (int iyy = std::max(0., iy - cutoff); iyy < iy + cutoff+1 && iyy < gridSize; iyy++) {
					for (int izz = std::max(0., iz - cutoff); izz < iz + cutoff+1 && izz < gridSize; izz++) {
						_map[((scrType*gridSize + izz)*gridSize + iyy )* gridSize + ixx] += 255 * exp256(-(
							(ixx * voxWidth - lx) * (ixx * voxWidth - lx) +
							(iyy * voxWidth - ly) * (iyy * voxWidth - ly) +
							(izz * voxWidth - lz) * (izz * voxWidth - lz)
							) / atomSpread2);
						
					}
				}
			}
			/////////////////////////////////////////////


			// R or L
			/////////////////////////////////////////////
			int rlType = 171;
			std::map<_pair,char>::iterator it_rl = _resList.find({atom->resSeqNumber,atom->chainId});
			if(it_rl != _resList.end())
			{
				switch (it_rl->second)
				{
				case 'R':
					rlType = 171;
					break;
				case 'L':
					rlType = 172;
					break;
				default:
					break;
				}
			}
			else{
				// std::cout<<"Notfound"<<std::endl;	
			}
			std::cout<<atom->resSeqNumber<<'\t'<<atom->chainId<<'\t'<<"rlType: "<<rlType<<'\t'<<"scrType: "<<scrType<<std::endl;
			
			for (int ixx = std::max(0., ix - cutoff); ixx < ix + cutoff+1 && ixx < gridSize; ixx++) {
				for (int iyy = std::max(0., iy - cutoff); iyy < iy + cutoff+1 && iyy < gridSize; iyy++) {
					for (int izz = std::max(0., iz - cutoff); izz < iz + cutoff+1 && izz < gridSize; izz++) {
						_map[((rlType*gridSize + izz)*gridSize + iyy )* gridSize + ixx] += 255 * exp256(-(
							(ixx * voxWidth - lx) * (ixx * voxWidth - lx) +
							(iyy * voxWidth - ly) * (iyy * voxWidth - ly) +
							(izz * voxWidth - lz) * (izz * voxWidth - lz)
							) / atomSpread2);
						
					}
				}
			}
			/////////////////////////////////////////////










			//std::cout << val << std::endl;
			return 1;
		}

	}
	return 0;
}


cVector3 cAtomicResidueMapper::addGradAtom(cAtom* atom, float* grads, float* typeVect) {
	cVector3 gradient(0, 0, 0);
	if (atom->elementSymbol[0] != 'C' && atom->elementSymbol[0] != 'O'&&  atom->elementSymbol[0] != 'N' && atom->elementSymbol[0] != 'S') {
		if (atom->elementSymbol[0] != 'H') {
			//std::cout << "Ignored atom " << atom << std::endl;
		}
		return gradient;
	}
	if (atom->residue->type() != cResidue::eType::AMINO) {
		return gradient;
	}
	int type = getType(atom);
	if (type >= 0) {
		cVector3 localPos = worldToLocal * (atom->getPosition() - origin);
		float gridLength = gridSize * voxWidth;
		int cutoff = 3;
		float cutoffDist = cutoff*atomSpread;
		if (localPos[0] > -cutoffDist && localPos[0] < gridLength + cutoffDist &&
			localPos[1] > -cutoffDist && localPos[1]  < gridLength + cutoffDist &&
			localPos[2] > -cutoffDist && localPos[2] < gridLength + cutoffDist) {
			int ix = localPos[0] / voxWidth;
			int iy = localPos[1] / voxWidth;
			int iz = localPos[2] / voxWidth;
			for (int ixx = std::max(0, ix - cutoff); ixx < ix + cutoff + 1 && ixx < gridSize; ixx++) {
				float dx = (ixx * voxWidth - localPos[0]);
				for (int iyy = std::max(0, iy - cutoff); iyy < iy + cutoff + 1 && iyy < gridSize; iyy++) {
					float dy = (iyy * voxWidth - localPos[1]);
					for (int izz = std::max(0, iz - cutoff); izz < iz + cutoff + 1 && izz < gridSize; izz++) {
						float dz = (izz * voxWidth - localPos[2]);
						float gradIntensity = 0;
						for (int k = 0; k < numRetype; k++) {
							gradIntensity += grads[((izz*gridSize + iyy)* gridSize + ixx) *numRetype + k] * typeVect[k];
						}
						gradient += -255*gradIntensity* exp(-(
							dx*dx +
							dy*dy +
							dz*dz
							) / atomSpread) * (dx*vX+dy*vY+dz*vZ);

					}
				}
			}
			//std::cout << val << std::endl;
			return gradient;
		}

	}
	return gradient;
}


int cAtomicResidueMapper::getType(cAtom* atom) {
	if (atom->residue != NULL && atom->residue->type() == cAminoResidue::eType::AMINO) {
		if (atomTypes.find(((cAminoResidue*) atom->residue)->aminoType()) != atomTypes.end()) {
			std::map<std::string,int> atomMap = atomTypes.at(((cAminoResidue*)atom->residue)->aminoType());
			if (atomMap.find(atom->name) != atomMap.end()) {
				return atomMap.at(atom->name);
			}
			else {
				if (atomTypes.at(cAminoResidue::AMINO_TYPE_END).find(atom->name) != atomTypes.at(cAminoResidue::AMINO_TYPE_END).end()) {
					return atomTypes.at(cAminoResidue::AMINO_TYPE_END).at(atom->name);
				}
				//std::cout << "Atom "<< atom->name <<" not found in the residue " << ((cAminoResidue*)atom->residue)->aminoType() << std::endl;
			}
		}
		else {
			std::cout << "Residue type not found " << ((cAminoResidue*) atom->residue)->aminoType() << std::endl;
		}
	}
	return -1;
}


int32_t* cProteinMapper::getMeta() {
	return _meta;
}

const uint8_t* cProteinMapper::getMap() {
	return _map;
}

void cResidueMapper::next() {
	++currentResIt;
}


int cCoarseResidueMapper::getType(cResidue* res) {
	if (res != NULL && res->type() == cAminoResidue::eType::AMINO) {
		cAminoResidue* ares = (cAminoResidue*)res;
		int type = ares->aminoType();
		if (type > 0 && type < residueTypes.size()) {
			return residueTypes[type];
		}
	}
	return -1;
}

int cCoarseResidueMapper::initResidueTypes() {
	// ALA
	residueTypes.push_back(0);
	// ARG
	residueTypes.push_back(1);
	// ASN
	residueTypes.push_back(2);
	// ASP
	residueTypes.push_back(3);
	// CYS
	residueTypes.push_back(4);
	// GLN
	residueTypes.push_back(5);
	// GLU
	residueTypes.push_back(6);
	// GLY
	residueTypes.push_back(7);
	// HIS
	residueTypes.push_back(8);
	// ILE
	residueTypes.push_back(9);
	// LEU
	residueTypes.push_back(10);
	// LYS
	residueTypes.push_back(11);
	// MET
	residueTypes.push_back(12);
	// PHE
	residueTypes.push_back(13);
	// PRO
	residueTypes.push_back(14);
	// SER
	residueTypes.push_back(15);
	// THR
	residueTypes.push_back(16);
	// TRP
	residueTypes.push_back(17);
	// TYR
	residueTypes.push_back(18);
	// VAL
	residueTypes.push_back(19);
	// SEC
	residueTypes.push_back(15);
	// MSE
	residueTypes.push_back(12);
	// AMINO_TYPE_END
	residueTypes.push_back(20);
	return 21;
}

cCoarseResidueMapper::cCoarseResidueMapper(size_t gs, float as, float vw,  bool orient) : cResidueMapper(gs, vw, orient), residueSpread(as)
{

	nbType = initResidueTypes();
	int metaSize = crmMetaSize();
	_meta = new int32_t[metaSize / sizeof(int32_t)];
	memset(_meta, 0, metaSize);
	_map = new uint8_t[gs*gs*gs*nbType*6];
	int headerSize = 8;
	header = new int[headerSize];
	header[0] = headerSize * sizeof(int);
	header[1] = 7919;
	header[2] = nbType*6;
	header[3] = gridSize;
	header[4] = gridSize;
	header[5] = gridSize;
	header[6] = metaSize;
	header[7] = -1;

}

cCoarseResidueMapper::~cCoarseResidueMapper() {

}

int cCoarseResidueMapper::addRes(cResidue* residue) {

	if (residue->type() != cResidue::eType::AMINO) {
		return 0;
	}
	cAminoResidue* ares = (cAminoResidue*)residue;
	int type = getType(ares);
	const cAtom*  cAlpha = ares->getAtomByLabel(cAminoResidue::eAtomLabel::CA);
	const cAtom*  carbon = ares->getAtomByLabel(cAminoResidue::eAtomLabel::C);
	const cAtom*  cBeta = ares->getAtomByLabel(cAminoResidue::eAtomLabel::CB);
	if (cAlpha == NULL) return 0;
	if (carbon == NULL) return 0;
	if (cBeta == NULL) {
		cBeta = ares->getAtomByLabel(cAminoResidue::eAtomLabel::HA1);
	}
	if (cBeta == NULL) return 0;
	cVector3 p1 = carbon->getPosition() - cAlpha->getPosition();
	cVector3 p2 = cBeta->getPosition() - cAlpha->getPosition();
	p1 = worldToLocal * p1;
	p2 = worldToLocal * p2;
	p1.normalize();
	p2.normalize();
	if (type >= 0) {
		cVector3 localPos = (cAlpha->getPosition() - origin);
		double lx = (localPos | vX), ly = (localPos | vY), lz = (localPos | vZ);
		float gridLength = gridSize * voxWidth;
		int cutoff = 3;
		float cutoffDist = cutoff*resSpread;
		if (lx > -cutoffDist && lx < gridLength + cutoffDist &&
			ly > -cutoffDist && ly < gridLength + cutoffDist &&
			lz > -cutoffDist && lz < gridLength + cutoffDist) {
			int ix = lx / voxWidth;
			int iy = ly / voxWidth;
			int iz = lz / voxWidth;
			for (int dirtype = 6*type; dirtype < 6*type+6; dirtype++) {
				float dirMul = (dirtype < 3 )? p1[dirtype] : p2[dirtype - 3];
				for (int ixx = std::max(0, ix - cutoff); ixx < ix + cutoff + 1 && ixx < gridSize; ixx++) {
					for (int iyy = std::max(0, iy - cutoff); iyy < iy + cutoff + 1 && iyy < gridSize; iyy++) {
						for (int izz = std::max(0, iz - cutoff); izz < iz + cutoff + 1 && izz < gridSize; izz++) {
							_map[((type*gridSize + izz)*gridSize + iyy)* gridSize + ixx] += 255 * dirMul * exp(-(
								(ixx * voxWidth - lx) * (ixx * voxWidth - lx) +
								(iyy * voxWidth - ly) * (iyy * voxWidth - ly) +
								(izz * voxWidth - lz) * (izz * voxWidth - lz)
								) / resSpread);

						}
					}
				}
			}
			//std::cout << val << std::endl;
			return 1;
		}

	}
	return 0;
}

void cCoarseResidueMapper::runCurrent() {
	memset(_map, 0, gridSize*gridSize *gridSize *nbType);
	setLocalFrame(&(*currentResIt), NULL, _orient);
	if (currentResIt->type() == cResidue::eType::AMINO) {
		cAminoResidue * aminoRes = (cAminoResidue*) &(*currentResIt);
		addRes(aminoRes);
		if (neighbourRes.find(aminoRes) != neighbourRes.end()) {
			std::vector<cResidue*> neighbour_residues = neighbourRes.at(aminoRes);
			int nbNeighb = 0;
			 #pragma omp parallel for
			for (const auto &residue_neighbours : neighbour_residues) {
				
				addRes(residue_neighbours);
			}
		}
		//std::cout << aminoRes->seqNumber << std::endl;
		((int32_t*)_meta)[0] = aminoRes->seqNumber;
		_meta[1] = aminoRes->getResCode();
		_meta[2] = -1000000;
	}
}
