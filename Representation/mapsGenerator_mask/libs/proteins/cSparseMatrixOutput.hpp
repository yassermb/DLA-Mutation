#pragma once

#include <string>
#include <vector>


// life cycle : initialize -> writeTriplet -> deinitialize
class cSparseMatrixOutput {
 public:
  virtual ~cSparseMatrixOutput() {}
  virtual bool initialize(const std::string &filename) = 0;
  virtual bool writeTriplet(size_t row, size_t col, double value) = 0;
  virtual bool deinitialize() = 0;
};

class cSparseMatrixOutputASCII : public cSparseMatrixOutput {
 public:
  cSparseMatrixOutputASCII() : file(NULL) {}
  ~cSparseMatrixOutputASCII();
  bool initialize(const std::string &filename);
  bool writeTriplet(size_t row, size_t col, double value);
  bool deinitialize();

 private:
  FILE *file;
};

class cSparseMatrixOutputMAT : public cSparseMatrixOutput {
 public:
  cSparseMatrixOutputMAT(size_t numRows, size_t numCols);
  ~cSparseMatrixOutputMAT();
  bool initialize(const std::string &filename);
  // FYI: writeTriplet column by column
  bool writeTriplet(size_t row, size_t col, double value);
  bool deinitialize();

 private:
  size_t numRows, numCols;
  std::vector<double> data;
  std::vector<int> ir;
  std::vector<int> jc;
  std::string filename;
};

bool getVectorFromFile(const std::string &file, std::vector<double> *vector);
