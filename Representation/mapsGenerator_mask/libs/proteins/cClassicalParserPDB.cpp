/*************************************************************************\

 Sergei Grudinin, 2010
 Mikhail Karasikov, 2014-2017
 All Rights Reserved.

 \**************************************************************************/

#include "cClassicalParserPDB.hpp"
#include "cAminoResidue.hpp"

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <fstream>

#ifdef TIMINGS
#include <ctime>
#endif


cClassicalParserPDB::eRecordType cClassicalParserPDB::readPDBLine(const std::string &line,
                                                                  cProtein *protein) {
  if (line.size() < 6)
    return eRecordType::PDB_EOF;

  if (!strncmp(line.data(), "ATOM ",  5)) {
    return eRecordType::PDB_ATOM;
  }
  else if (!strncmp(line.data(), "HETATM", 6)) {
    return eRecordType::PDB_HETATM;
  }
  else if (!strncmp(line.data(), "TER", 3)) {
    return eRecordType::PDB_TER;
  }
  else if (!strncmp(line.data(), "ENDMDL", 6)) {
    return eRecordType::PDB_ENDMDL;
  }
  else if (!strncmp(line.data(), "END", 3)) {
    return eRecordType::PDB_END;
  }
  else if (!strncmp(line.data(), "REMARK", 6)) {
    char tmpStr28[28];
    tmpStr28[27]='\0';
    strncpy(tmpStr28, line.data(), 27);
    if (!strncmp(tmpStr28,"REMARK 465   M RES C SSSEQI",27)) {
      switch_PDB_MISSING_RESIDUES = true;
#if DEBUG
      std::cout<<"Missing Residues:\n";
      std::cout<<tmpStr28 + 10<<std::endl;
#endif
    } else if (!strncmp(tmpStr28,"REMARK 465",10) && switch_PDB_MISSING_RESIDUES) {
#if DEBUG
      std::cout<<tmpStr28 + 10<<std::endl;
#endif
      return eRecordType::PDB_MISSING_RESIDUES;
    } else if (!strncmp(line.data(), "REMARK 500  M RES CSSEQI        PSI       PHI", 45)) {
#if DEBUG
      switch_PDB_EXPECTED_PHI_PSI = true;
#endif
    } else if (!strncmp(line.data(), "REMARK 500", 10) && switch_PDB_EXPECTED_PHI_PSI) {
      if (line.size() < 45 || line[18] == ' ') {
        switch_PDB_EXPECTED_PHI_PSI = false;
      } else {
        size_t resSeqNumber;
        double phi, psi;
        if (sscanf(line.data() + 19, "%zd%lf%lf", &resSeqNumber, &psi, &phi) != 3) {
          std::cerr << "Error reading expected PHI, PSI:\n" << line << std::endl;
          exit(1);
        }
        protein->vExpectedPhiPsi.push_back({line[18], resSeqNumber, phi, psi});
      }
    }
    return eRecordType::PDB_REMARK;
  }
  return eRecordType::PDB_UNKNOWN;
}

cAtom* cClassicalParserPDB::parsePDBAtom(const std::string &field) {
  cAtom *atom = new cAtom();
  if (atom->parseFromPDBField(field)) {
    return atom;
  }
  std::cerr << "Error when parsing Atom from PDB field:\n" << field << std::endl;
  delete atom;
  return NULL;
}

void writeLinePDB(FILE *pdbFile, const cAtom &atom) {
  atom.toFile(pdbFile);
}

bool cClassicalParserPDB::parse(const std::string &pdbName,
                                cProtein *protein,
                                std::function<bool(const cAtom &)> skip_atom,
                                bool skip_structure_errors,
                                bool skip_pdb_errors) {
  auto &chains = protein->chains_;
  for (const auto *chain : chains) {
    delete chain;
  }
  chains.clear();

  protein->vExpectedPhiPsi.clear();

  // ******************
  // parse the pdb file
  // ******************

  std::ifstream ifs(pdbName);
  std::string line;
  if (!ifs.good()) {
    std::cerr << "Error: could not open pdb file " << pdbName.c_str() << "\n";
    return false;
  }

#if DEBUG
  std::cout << "Reading PDB file " << pdbName.c_str() << "...\n";
#endif
#ifdef TIMINGS
  std::clock_t tp0 = std::clock();
#endif

  cProteinChain *curChain = NULL;
  cResidue *curResidue = NULL;
  cResidue *prevResidue = NULL;

  // read the next line
  while (std::getline(ifs, line)) {
    eRecordType indx = readPDBLine(line, protein);

    if (indx == PDB_END && indx == PDB_EOF)
      break;

    // if (indx == PDB_TER) {
      // // new chain
      // curChain = NULL;
    // }

    //if (indx == PDB_ATOM || indx == PDB_HETATM) {
    if (indx != PDB_ATOM)
      continue;

    cAtom *curAtom = parsePDBAtom(line);
    if (!curAtom) {
      if (skip_pdb_errors) {
        continue;
      } else {
        return false;
      }
    }

    if (skip_atom(*curAtom)
        || curAtom->insertion != ' '
        || !strcmp(curAtom->name, "CEN")) {
      delete curAtom;
      continue;
    }

    if (!curChain || curAtom->chainId != curChain->chainId) {
		if (curChain && curResidue) {
			addResidueIfDefined(curResidue, skip_structure_errors, curChain);

		}
      // new chain
      curChain = new cProteinChain(curAtom->chainId);
	  curResidue = NULL;
      chains.push_back(curChain);
    }
    if (!curResidue || curResidue->chainId != curAtom->chainId ||
                       curResidue->seqNumber != curAtom->resSeqNumber) {
      // next residue
      addResidueIfDefined(curResidue, skip_structure_errors, curChain);
      curResidue = cResidue::create(curAtom->resName, curAtom->chainId,
                                    curAtom->segId, curAtom->resSeqNumber, curResidue);
    }
    try {
      curResidue->addAtom(curAtom);
    } catch (std::exception &e) {
      std::cerr << "Error: " << e.what() << std::endl;
      if (skip_structure_errors) {
        continue;
      } else {
        return false;
      }
    } catch (...) {
      std::cerr << "Unknown exception" << std::endl;
      if (skip_structure_errors) {
        continue;
      } else {
        return false;
      }
    }
  }
  assert(!curResidue || curChain);
  addResidueIfDefined(curResidue, skip_structure_errors, curChain);

#ifdef TIMINGS
  std::cout << "File is read in " << (std::clock() - tp0) / (double)CLOCKS_PER_SEC
            << " seconds\n";
#endif

  if (!chains.size()) {
    std::cerr << "Error: No chains have been read\n";
    return false;
  }

  if (!protein->numAtoms()) {
    std::cerr << "Error: No atoms in the protein\n";
    return false;
  }

#if DEBUG
  // std::cout << "OK (" << atoms.size() << " atoms) \n";
#endif
  return true;
}

bool cClassicalParserPDB::addResidueIfDefined(cResidue *residue,
                                              bool skip_structure_errors,
                                              cProteinChain *chain) {
  if (!residue)
    return false;

  auto *res = dynamic_cast<const cAminoResidue *>(residue);

  if (res && !res->isBackboneValid() && skip_structure_errors) {
    delete residue;
    return false;
  }
  chain->addResidue(residue);
  return true;
}
