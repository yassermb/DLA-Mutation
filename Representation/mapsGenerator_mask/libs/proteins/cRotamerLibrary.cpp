#include "cRotamerLibrary.hpp"

#include <assert.h>
#include <fstream>

const size_t kAminoTypesNumber = 20;  // Number of types of amino acids
const size_t kBBAnNum = 37;  // Number of different angles phi in library

cPDBDunbrackBDLib::cPDBDunbrackBDLib()
  : isInitialized(false), libData(kAminoTypesNumber * kBBAnNum * kBBAnNum + 1) {}

bool cPDBDunbrackBDLib::readStandartLib(const std::string &libName) {
  std::ifstream ifs(libName);
  std::string line;
  if (!ifs.good())
    return false;
  while (std::getline(ifs, line)) {
    if (line[0] == '#' || line.size() < 78)
      continue;

    std::map<std::string, cAminoResidue::eAminoType>::const_iterator it
      = cAminoResidue::mAminoType.find(line.substr(0,3));
    if (it == cAminoResidue::mAminoType.end())
      continue;

    size_t curAmino = static_cast<size_t>(it->second);

    int phi, psi, r[4];
    cRotamer rotamer;

    if (sscanf((line+='\0').data() + 3, "%d%d%*d%d%d%d%d%lf%lf%lf%lf%lf",
               &phi, &psi, &r[0], &r[1], &r[2], &r[3], &rotamer.prob,
               &rotamer.chi[0], &rotamer.chi[1], &rotamer.chi[2], &rotamer.chi[3] ) < 11)
      return false;

    int curPhi = phi/10 + 18,
        curPsi = psi/10 + 18;
    libData[curAmino * kBBAnNum * kBBAnNum + curPhi * kBBAnNum + curPsi].push_back(rotamer);
  }
  libData[kAminoTypesNumber * kBBAnNum * kBBAnNum].push_back({{0, 0, 0, 0}, 1.0});
  return true;
}

bool cPDBDunbrackBDLib::writeFormattedLib(const std::string &libName) {
  std::ofstream ofs(libName);
  if (!ofs.good())
    return false;

  const int kBufferSize = 100;
  char buffer[kBufferSize];
  for (size_t i = 0; i < libData.size(); ++i) {
    if (!libData[i].size())
      continue;

    if (snprintf(buffer, kBufferSize, "%-6zd%-3zd", i, libData[i].size()) != 9)
      return false;

    ofs << buffer;
    for (size_t j = 0; j < libData[i].size(); ++j) {
      if (snprintf(buffer, kBufferSize, "%-6.1lf%-6.1lf%-6.1lf%-6.1lf%-8.6lf",
                              libData[i][j].chi[0], libData[i][j].chi[1],
                              libData[i][j].chi[2], libData[i][j].chi[3],
                              libData[i][j].prob) != 32)
        return false;

      if (!(ofs << buffer))
        return false;
    }
    ofs << "\n";
  }
  return ofs.good();
}

bool cPDBDunbrackBDLib::readFormattedLib(const std::string &libName) {
  std::ifstream ifs(libName);
  if (!ifs.good())
    return false;

  size_t i, n;
  double chi[4], prob;
  std::string line;
  while (std::getline(ifs, line)) {
    if (sscanf(line.data(), "%6zd%3zd", &i, &n) < 2)
      return false;

    const char *pos = line.data() + 9;
    for (size_t j = 0; j < n; ++j) {
      if (sscanf(pos, "%6lf%6lf%6lf%6lf%8lf",
                              &chi[0], &chi[1], &chi[2], &chi[3], &prob) < 5)
        return false;
      libData[i].push_back({{chi[0], chi[1], chi[2], chi[3]}, prob});
      pos += 32;
    }
  }
  return ifs.eof();
}

bool cPDBDunbrackBDLib::writeBinaryLib(const std::string &libName) {
  std::ofstream ofs(libName, std::ofstream::binary);
  if (!ofs.good())
    return false;

  size_t libDataSize = libData.size();
  if (!ofs.write(reinterpret_cast<const char *>(&libDataSize), sizeof(libDataSize)))
    return false;

  for (size_t i = 0; i < libDataSize; ++i) {
    size_t size = libData[i].size();
    if (!ofs.write(reinterpret_cast<const char *>(&size), sizeof(size)))
      return false;
    if (!ofs.write(reinterpret_cast<const char *>(libData[i].data()),
                   sizeof(cRotamer) * size))
      return false;
  }
  return ofs.good();
}

bool cPDBDunbrackBDLib::readBinaryLib(const std::string &libName) {
  std::ifstream ifs(libName, std::ofstream::binary);
  if (!ifs.good())
    return false;

  size_t libDataSize;
  if (!ifs.read(reinterpret_cast<char *>(&libDataSize), sizeof(libDataSize)))
    return false;

  libData.resize(libDataSize);
  for (size_t i = 0; i < libDataSize; ++i) {
    size_t size;
    if (!ifs.read(reinterpret_cast<char *>(&size), sizeof(size)))
      return false;
    libData[i].resize(size);
    if (!ifs.read(reinterpret_cast<char *>(libData[i].data()),
                  sizeof(cRotamer) * size))
      return false;
  }
  return ifs.good();
}

bool cPDBDunbrackBDLib::initialize(const std::string &libName) {
  if (readBinaryLib(libName + ".bin"))
    return isInitialized = true;

  remove((libName + ".bin").c_str());
  if (readFormattedLib(libName + ".formatted")) {
    writeBinaryLib(libName + ".bin");
    return isInitialized = true;
  }
  remove((libName + ".formatted").c_str());
  if (readStandartLib(libName)) {
    writeFormattedLib(libName + ".formatted");
    writeBinaryLib(libName + ".bin");
    return isInitialized = true;
  }
  return false;
}

const std::vector<cRotamer> *
        cPDBDunbrackBDLib::getRotamers(const cAminoResidue &residue) const {
  assert(isInitialized);
  double phi = residue.getPhi(), psi = residue.getPsi();
  if (phi == 181.0)
    phi = -60.0;

  if (psi == 181.0)
    psi = 60.0;

  assert(phi >= -180.1 && phi <= 180.1 && psi >= -180.1 && psi <= 180.1);

  if (residue.aminoType() == cAminoResidue::eAminoType::ALA
      || residue.aminoType() == cAminoResidue::eAminoType::GLY)
    return &libData[kAminoTypesNumber * kBBAnNum * kBBAnNum];

  size_t curAmino = static_cast<size_t>(residue.aminoType());
  int curPhi = static_cast<int>(round(phi/10)) + 18,
      curPsi = static_cast<int>(round(psi/10)) + 18;
  return &libData[curAmino * kBBAnNum * kBBAnNum + curPhi * kBBAnNum + curPsi];
}
