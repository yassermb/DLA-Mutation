/*************************************************************************\

Sergei Grudinin, 2010
Mikhail Karasikov, 2014-2017
All Rights Reserved.

\**************************************************************************/

#pragma once

#include "cProtein.hpp"


class cClassicalParserPDB {
 public:
	cClassicalParserPDB()
    : switch_PDB_MISSING_RESIDUES(false), switch_PDB_EXPECTED_PHI_PSI(false) {};

  // Missing of records TER is allowed
	bool parse(const std::string &pdbName,
             cProtein *protein,
             std::function<bool(const cAtom &)> skip_atom = [](const cAtom &) { return false; },
             bool skip_structure_errors = false,
             bool skip_pdb_errors = false);
 private:
	//  record types
	enum eRecordType {
		PDB_HEADER, PDB_REMARK, PDB_ATOM, PDB_HETATM, PDB_CONECT, PDB_UNKNOWN, PDB_END,
    PDB_EOF, PDB_CRYST1, BIOMT, PDB_TER, SEQADV, SEQRES, PDB_MISSING_RESIDUES, PDB_ENDMDL
	};

	//enum eFileType {PDB, PQR, XYZ, NONE};

	bool    switch_PDB_MISSING_RESIDUES;
	bool    switch_PDB_EXPECTED_PHI_PSI;

  eRecordType readPDBLine(const std::string &line, cProtein *protein);
	static cAtom* parsePDBAtom(const std::string &field);
  static bool addResidueIfDefined(cResidue *residue,
                                  bool skip_structure_errors,
                                  cProteinChain *chain);
};

void writeLinePDB(FILE* pdbFile, const cAtom &curAtom);
