#pragma once

#include "cAminoResidue.hpp"


class cRotamerLibrary {
 public:
  virtual ~cRotamerLibrary() {}
  virtual bool initialize(const std::string &libName) = 0;
  // vector of the most probable rotamers sorted by the probability
  virtual const std::vector<cRotamer> *
                  getRotamers(const cAminoResidue &residue) const = 0;
};

// Dunbrack back-bone dependent rotamer library
class cPDBDunbrackBDLib : public cRotamerLibrary {
 public:
  cPDBDunbrackBDLib();
  // ALL.bbdep.rotamers.lib
  bool initialize(const std::string &libName);

  const std::vector<cRotamer> *
          getRotamers(const cAminoResidue &residue) const;

 private:
  bool isInitialized;
  std::vector<std::vector<cRotamer>> libData;
  bool readStandartLib(const std::string &libname);
  bool writeFormattedLib(const std::string &libname);
  bool readBinaryLib(const std::string &libname);
  bool writeBinaryLib(const std::string &libname);
  bool readFormattedLib(const std::string &libname);
};
