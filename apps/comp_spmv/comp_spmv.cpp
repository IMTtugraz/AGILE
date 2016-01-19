// Copyright (C) 2010-2011 Institute of Medical Engineering,
// Graz University of Technology
//
// This program is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
// version.
//
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License along with
// this program; if not, see <http://www.gnu.org/licenses>.

// $Id: comp_spmv.cpp 476 2011-06-16 08:54:14Z freiberger $

#include "agile/gpu_environment.hpp"
#include "agile/network/gpu_communicator.hpp"
#include "agile/gpu_vector.hpp"
#include "agile/gpu_cs_matrix.hpp"

#include "../../test/test_defines.h"

#include "./spmv/SpMV.h"
#include "./config.h"

#include <stdio.h>
#include <time.h>
#include <sys/time.h>

bool fileIsReadable(const char *filename)
{
  FILE* f;
  if ((f = fopen(filename, "r")) == NULL) {
    return false;
  }
  fclose(f);
  return true;
}

void readSparseMatrix(const char *filename, int format, unsigned &numRows,
  unsigned &numCols, std::vector<unsigned> &nnz, std::vector<unsigned> &index,
  std::vector<float> &data)
{
  FILE *f;
  f = fopen(filename,"r");
  if (!f) {
      fprintf(stderr,"Cannot open file: %s\n",filename);
      exit(-1);
  }

  char line[256];

  // skip matrix market file header
  while ( (fgets(line, 256, f)) != NULL) {
    if (line[0] != '%') break;
  }

  unsigned numNZEntries;
  if ( (sscanf(line,"%d %d %d", &numRows, &numCols, &numNZEntries)) != 3) {
    exit(-1);
  }

  NZEntry *nzEntries = (NZEntry *) malloc(sizeof(NZEntry) * numNZEntries);

  // resize vectors
  index.resize(numNZEntries, 0);
  data.resize(numNZEntries, 0);
  // format 1 .. column major, 0 .. row major
  nnz.resize(format == 1 ? numCols : numRows, 0);

  NZEntry e;
  for (unsigned i=0; i<numNZEntries; ++i)
  {
    if (fscanf(f,"%d %d %f\n", &(e.rowNum), &(e.colNum), &(e.val)) != 3) {
      exit(-1);
    }
    e.rowNum--; e.colNum--; // zero based indizes
    nzEntries[i] = e;
  }

  // finished reading from file, so we close it
  fclose(f);

  // sort into row-major order or column major order based on the format
  // and copy values to result vectors
  if (format == 0) // ROW-MAJOR
  {
    qsort(nzEntries, numNZEntries, sizeof(NZEntry), cmpRow);
    for (unsigned i=0; i<numNZEntries; ++i) {
      nnz[nzEntries[i].rowNum]++;
      index[i] = nzEntries[i].colNum;
      data[i] = nzEntries[i].val;
    }
  }
  else if (format == 1)  // COLUMN-MAJOR
  {
    qsort(nzEntries, numNZEntries, sizeof(NZEntry), cmpCol);
    for (unsigned i=0; i<numNZEntries; ++i) {
      nnz[nzEntries[i].colNum]++;
      index[i] = nzEntries[i].rowNum;
      data[i] = nzEntries[i].val;
    }
  }

  free(nzEntries);
}











int main(int argc, char** argv)
{
  // timer
  struct timeval st, et;
  float gputime = 0.0, cputime = 0.0;

  // read Sparse Matrix from file or generate
  if (argc < 2 || argc > 4) {
      printf("Correct Usage: <executable> <input matrix file>\n");
      exit(-1);
  }

  // init the network
  agile::NetworkEnvironment environment(argc, argv);

  // allocate a GPU
  typedef agile::GPUCommunicator<unsigned, float, float> communicator_type;
  communicator_type com;
  com.allocateGPU();

  char spmfileName[256];
  strcpy(spmfileName, argv[1]);
  if (!fileIsReadable(spmfileName))
  {
    printf("Non-existent input matrix file\n");
    exit(-1);
  }

  unsigned m_num_rows, m_num_cols;
  std::vector<unsigned> m_row_nnz;
  std::vector<unsigned> m_column_index;
  std::vector<float> m_data;

  // read in matrix from matrix-market file
  readSparseMatrix(spmfileName, 0, m_num_rows, m_num_cols, m_row_nnz,
    m_column_index, m_data);

  std::cout << m_num_rows << "\t" << m_num_cols << "\t";
/*
  PRINT_VEC("m_row_nnz", m_row_nnz);
  PRINT_VEC("m_column_index", m_column_index);
  PRINT_VEC("m_data", m_data);
*/

  // init gpu matrix
  agile::GPUCSMatrix<float> A(m_row_nnz, m_column_index, m_data);

  // init random vector
  std::vector<float> x_host(m_num_cols, 0);
  srand(time(NULL));
  for (unsigned i=0; i<m_num_cols; ++i)
    x_host[i] = rand() / (float)RAND_MAX;

//PRINT_VEC("RANDOM X VECTOR", x_host);

  // init gpu vector
  agile::GPUVector<float> x(m_num_cols);
  x.assignFromHost(x_host.begin(), x_host.end());

  // init result gpu vector: y
  agile::GPUVector<float> y(m_num_rows);

  // start time
  gettimeofday(&st, NULL);

  for (unsigned t=0; t<NUM_ITER; ++t)
  {
    // gpu multiplication
    agile::multiply(A, x, y);

    cudaThreadSynchronize();
  }

  // stop time
  gettimeofday(&et, NULL);
  gputime = ((et.tv_sec-st.tv_sec)*1000.0 + (et.tv_usec - st.tv_usec)/1000.0)/NUM_ITER;

  // transfer GPU multiplication result back to cpu
  std::vector<float> y_host;
  y.copyToHost(y_host);


  //----------------- CPU computation from ibm demo ---------------------------
  SpMatrix m;
  readSparseMatrix(&m, spmfileName, 0);
  unsigned int numNonZeroElements = m.numNZEntries;
  unsigned int memSize_row = sizeof(float) * m_num_rows;

  // allocate host memory
  float* h_x = (float*) malloc(memSize_row); 

  #if PADDED_CSR
    float *h_val;
    unsigned int *h_indices, *h_rowIndices;
    genPaddedCSRFormat(&m, &h_val, &h_rowIndices, &h_indices);
  #else
    float* h_val = (float*) malloc(sizeof(float)*numNonZeroElements);
    unsigned int* h_indices = (unsigned int*) malloc(sizeof(int)*numNonZeroElements);
    unsigned int* h_rowIndices = (unsigned int*) malloc(sizeof(int)*(m_num_rows+1));
    genCSRFormat(&m, h_val, h_rowIndices, h_indices);
  #endif

  // CPU REFERENCE
  float* reference = (float*) malloc(memSize_row);
#if EXEC_CPU
  #if TIMER
  gettimeofday(&st, NULL);
  #endif
  // compute reference solution
  #if BCSR
  float *val;
  unsigned int *rowIndices, *indices;
  unsigned int numblocks;
  genBCSRFormat(&m, &val, &rowIndices, &indices, &numblocks, BCSR_r, BCSR_c);
  computeSpMV_BCSR(reference, val, rowIndices, indices, &(x_host[0]), m_num_rows, m_num_cols, BCSR_r, BCSR_c);
  #else
  computeSpMV(reference, h_val, h_rowIndices, h_indices, &(x_host[0]), m_num_rows);
  #endif
  #if TIMER
  gettimeofday(&et, NULL);
  cputime = (et.tv_sec-st.tv_sec)*1000.0 + (et.tv_usec - st.tv_usec)/1000.0;
  #endif
#endif

  float flops= ((numNonZeroElements * 2) / (gputime*1000000));
  //printf("GPU (ms) \tCPU (ms) \tGFLOPS\n");
  printf("%f\t%f\t%f\t", gputime, cputime, flops);

#if VERIFY
  // check result
  float error_norm, ref_norm, diff;
  error_norm = 0;
  ref_norm = 0;
  for (unsigned i = 0; i < m_num_rows; ++i) {
      diff = reference[i] - y_host[i];
      error_norm += diff * diff;
      ref_norm += reference[i] * reference[i];
  }
  error_norm = (float)sqrt((double)error_norm);
  ref_norm = (float)sqrt((double)ref_norm);

  if (fabs(ref_norm) < 1e-7)
    printf ("Test FAILED");
  else
    printf( "Test %s", ((error_norm / ref_norm) < 1e-6f) ? "PASSED" : "FAILED");

#endif

  free(reference);
  free(h_x);
  #if !PADDED_CSR
    free(h_val);
    free(h_indices);
    free(h_rowIndices);
  #endif

  return 0;
}





// End of $Id: comp_spmv.cpp 476 2011-06-16 08:54:14Z freiberger $.

