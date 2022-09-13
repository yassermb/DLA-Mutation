#include "cSparseMatrixOutput.hpp"

#include <assert.h>
#include <stdio.h>

#define NO_MATIO 1

#if !NO_MATIO
#ifndef _WIN32
#include <matio.h>

#endif // !_WIN32
#endif // !NO_MATIO



cSparseMatrixOutputASCII::~cSparseMatrixOutputASCII() {
  if (file)
    deinitialize();
}

bool cSparseMatrixOutputASCII::initialize(const std::string &filename) {
  assert(!file);
  if ((file = fopen((filename + ".txt").c_str(), "w")))
    return true;

  perror("cSparseMatrixOutputASCII error");
  return false;
}

bool cSparseMatrixOutputASCII::writeTriplet(size_t row, size_t col, double value) {
  assert(file);
  if (fprintf(file, "%zd %zd %lf\n", row + 1, col + 1, value) > 0)
    return true;

  perror("cSparseMatrixOutputASCII error");
  return false;
}

bool cSparseMatrixOutputASCII::deinitialize() {
  assert(file);
  bool ret = !fclose(file);
  file = NULL;
  return ret;
}


cSparseMatrixOutputMAT::cSparseMatrixOutputMAT(size_t numRows_, size_t numCols_)
    : numRows(numRows_), numCols(numCols_) {
  assert(numRows > 0 && numCols > 0);
}

cSparseMatrixOutputMAT::~cSparseMatrixOutputMAT() {
  if (filename.size())
    deinitialize();
}

bool cSparseMatrixOutputMAT::initialize(const std::string &filename_) {
  filename = filename_;
  return filename.size();
}

bool cSparseMatrixOutputMAT::writeTriplet(size_t row, size_t col, double value) {
  assert(row < numRows && col < numCols && col + 1 >= jc.size());
  while (col != jc.size() - 1)
    jc.push_back(data.size());

  ir.push_back(row);
  data.push_back(value);
  return true;
}

bool cSparseMatrixOutputMAT::deinitialize() {
	#if !NO_MATIO
#ifndef _WIN32
  assert(filename.size());

  while (jc.size() != numCols + 1)
    jc.push_back(data.size());

  mat_t *mat = Mat_CreateVer((filename + ".mat").c_str(), NULL, MAT_FT_DEFAULT);
  if(!mat) {
    fprintf(stderr, "ERROR: cSparseMatrixOutputMAT cant create MAT file\n");
    return false;
  }

  mat_sparse_t sparse_matrix = { static_cast<int>(data.size()),
                                 ir.data(),
                                 static_cast<int>(ir.size()),
                                 jc.data(),
                                 static_cast<int>(jc.size()),
                                 static_cast<int>(data.size()),
                                 data.data() };
  size_t dims[2] = { numRows, numCols };

  std::string varName = "SpMat_";
  varName += filename.substr(filename.find_last_of(".\\/") + 1);
  matvar_t *matvar = Mat_VarCreate(varName.c_str(), MAT_C_SPARSE, MAT_T_DOUBLE,
                                   2, dims, (void *)&sparse_matrix, 0);
  Mat_VarWrite(mat, matvar, MAT_COMPRESSION_ZLIB);

  Mat_VarFree(matvar);
  Mat_Close(mat);
  data.clear();
  ir.clear();
  jc.clear();
  filename.clear();
  return true;
#else
	fprintf(stderr, "ERROR: unsupported MAT file on WINDOWS \n");
	return false;
#endif 
#else
	fprintf(stderr, "ERROR: compiled without matio support \n");
	return false;

#endif
}

bool getVectorFromFile(const std::string &filename, std::vector<double> *vector) {
  vector->clear();
  FILE *file = fopen(filename.c_str(), "r");
  if (!file) {
    perror("getVectorFromFile error");
    return false;
  }
  double value;
  while (fscanf(file, "%lf", &value) == 1)
    vector->push_back(value);

  if (feof(file)) {
    if (fclose(file)) {
      perror("getVectorFromFile error");
      return false;
    }
    return true;
  }
  perror("getVectorFromFile error");
  fclose(file);
  return false;
}
