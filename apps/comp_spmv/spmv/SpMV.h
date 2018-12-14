/*
 * IBM Sparse Matrix-Vector Multiplication Toolkit for Graphics Processing Units
 * (c) Copyright IBM Corp. 2008, 2009.  All Rights Reserved.
 */ 

#ifndef __SPMV_H__
#define __SPMV_H__

#define ceild(n,d)  ceil(((float)(n))/((float)(d)))
#define floord(n,d) floor(((float)(n))/((float)(d)))
#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))

struct nzInfo
{
    unsigned int rowNum;
    unsigned int colNum;
    float val;
};

typedef struct nzInfo NZEntry;

struct SpM
{
    unsigned int numRows;
    unsigned int numCols;
    unsigned int numNZEntries;
    NZEntry *nzentries;
    unsigned int *rowPtrs;
    unsigned int *colPtrs;
};

typedef struct SpM SpMatrix;

struct SpMGPU
{
    float* d_val;
    unsigned int* d_indices;
    unsigned int* d_rowIndices;
    unsigned int* d_ins_indices;
    unsigned int* d_ins_rowIndices;
    unsigned int* d_ins_inputList;
};

typedef struct SpMGPU SpMatrixGPU;

int cmpRow(const void *e1, const void *e2);
int cmpCol(const void *e1, const void *e2);
void readInputVector(float *y, const char *filename, int numCols);
void writeOutputVector(float *x, const char *filename, int numRows);
void readSparseMatrix(SpMatrix *m, const char *filename, int format);
void genCSRFormat(SpMatrix *m, float *val, unsigned int *rowIndices, unsigned int *indices);
void genCSCFormat(SpMatrix *m, float *val, unsigned int *colIndices, unsigned int *indices);
void genBCSRFormat(SpMatrix *m, float **val, unsigned int **rowIndices, unsigned int **indices,
                   unsigned int *numblocks, unsigned int bsx, unsigned int bsy);
void genPaddedCSRFormat(SpMatrix *m, float **val, unsigned int **rowIndices, unsigned int **indices);
void allocateSparseMatrixGPU(SpMatrixGPU *spm, SpMatrix *m, float *h_val, unsigned int *h_rowIndices,
				unsigned int *h_indices, const unsigned int numRows, const unsigned int numCols);
void SpMV_cuda(float *x, SpMatrixGPU *spm, const float *y, const unsigned int numRows,
               const unsigned int numCols, const unsigned int numNonZeroElements);
void computeSpMV(float *x, const float *val, const unsigned int *rowIndices, const unsigned int *indices,
                 const float *y, const unsigned int numRows);
void computeSpMV_BCSR(float *x, const float *val, const unsigned int *rowIndices,
                 const unsigned int *indices, const float *y, const unsigned int numRows,
                 const unsigned int numCols, const unsigned int bsx, const unsigned int bsy);

#endif /* __SPMV_H__ */
