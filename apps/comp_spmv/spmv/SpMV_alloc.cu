/*
 * IBM Sparse Matrix-Vector Multiplication Toolkit for Graphics Processing Units
 * (c) Copyright IBM Corp. 2008, 2009.  All Rights Reserved.
 */ 

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include "../config.h"
#include "SpMV.h"
#include "SpMV_inspect.h"

void allocateSparseMatrixGPU(SpMatrixGPU *spm, SpMatrix *m, float *h_val, unsigned int *h_rowIndices,
               unsigned int *h_indices, const unsigned int numRows, const unsigned int numCols) 
{

    unsigned int numNonZeroElements = h_rowIndices[numRows];
    unsigned int memSize_val = sizeof(float) * numNonZeroElements;
    unsigned int memSize_row = sizeof(int) * numRows;
    //unsigned int memSize_col = sizeof(float) * numCols;

  #if INSPECT
    unsigned int *ins_rowIndices, *ins_indices, *h_indicesFill, *h_rowIndicesFill;
    unsigned int nnz_fill, ins_numblocks, *ins_nnzCount_block, *ins_yCount_block;
    float *h_valFill;
    int insBStat;

    #if VAR_BLOCK
    insBStat = inspectVarBlock(m, &h_valFill, &h_indicesFill, &h_rowIndicesFill, &ins_rowIndices, &ins_indices, &ins_numblocks, &ins_nnzCount_block, &ins_yCount_block, &nnz_fill, (BLOCKSIZE/HALFWARP), INSPECT_BLOCK_c, VAR_COLUMN);
    if (insBStat == ERR_INSUFFICIENT_MEMORY) { printf("Insufficient Memory while malloc in inspectVarBlock\n"); exit(-1); }
    #else
    insBStat = inspectBlock(m, &ins_rowIndices, &ins_indices, &ins_numblocks, &ins_nnzCount_block, &ins_yCount_block, (BLOCKSIZE/HALFWARP), INSPECT_BLOCK_c);
    if (insBStat == ERR_INSUFFICIENT_MEMORY) { printf("Insufficient Memory while malloc in inspectBlock\n"); exit(-1); }
    #endif
  #endif

  #if INSPECT_INPUT
    unsigned int *ins_rowIndices, *ins_indices, *ins_inputList; 
    unsigned int ins_numblocks, ins_inputListCount;
    inspectInputBlock(m, &ins_inputList, &ins_rowIndices, &ins_indices, &ins_numblocks, &ins_inputListCount, (BLOCKSIZE/HALFWARP),HALFWARP);
  #endif

  #if (INSPECT && VAR_BLOCK && C_GLOBAL_OPT)
    cudaMalloc((void**) &(spm->d_val), sizeof(float)*nnz_fill);
    cudaMalloc((void**) &(spm->d_indices), sizeof(int)*nnz_fill);
  #else
    cudaMalloc((void**) &(spm->d_val), memSize_val);
    cudaMalloc((void**) &(spm->d_indices), sizeof(int)*numNonZeroElements);
  #endif
    cudaMalloc((void**) &(spm->d_rowIndices), memSize_row+sizeof(int));
  #if INSPECT
    cudaMalloc((void**) &(spm->d_ins_indices), sizeof(int)*ins_numblocks);
    cudaMalloc((void**) &(spm->d_ins_rowIndices), sizeof(int)*(1+(int)ceild(numRows,(BLOCKSIZE/HALFWARP))));
  #endif
  #if INSPECT_INPUT
    cudaMalloc((void**) &(spm->d_ins_indices), sizeof(int)*(ins_numblocks+1));
    cudaMalloc((void**) &(spm->d_ins_rowIndices), sizeof(int)*(1+(int)ceild(numRows,(BLOCKSIZE/HALFWARP))));
    cudaMalloc((void**) &(spm->d_ins_inputList), sizeof(int)*ins_inputListCount);
  #endif

  #if (INSPECT && VAR_BLOCK && C_GLOBAL_OPT)
    cudaMemcpy(spm->d_val, h_valFill, sizeof(float)*nnz_fill, cudaMemcpyHostToDevice);
    cudaMemcpy(spm->d_indices, h_indicesFill, sizeof(int)*nnz_fill, cudaMemcpyHostToDevice);
    cudaMemcpy(spm->d_rowIndices, h_rowIndicesFill, memSize_row+sizeof(int), cudaMemcpyHostToDevice);
  #else
    cudaMemcpy(spm->d_val, h_val, memSize_val, cudaMemcpyHostToDevice);
    cudaMemcpy(spm->d_indices, h_indices, sizeof(int)*numNonZeroElements, cudaMemcpyHostToDevice);
    cudaMemcpy(spm->d_rowIndices, h_rowIndices, memSize_row+sizeof(int), cudaMemcpyHostToDevice);
  #endif
  #if INSPECT
    cudaMemcpy(spm->d_ins_indices, ins_indices, sizeof(int)*ins_numblocks, cudaMemcpyHostToDevice);
    cudaMemcpy(spm->d_ins_rowIndices, ins_rowIndices, sizeof(int)*(1+(int)ceild(numRows,(BLOCKSIZE/HALFWARP))), cudaMemcpyHostToDevice);
  #endif
  #if INSPECT_INPUT
    cudaMemcpy(spm->d_ins_indices, ins_indices, sizeof(int)*(ins_numblocks+1), cudaMemcpyHostToDevice);
    cudaMemcpy(spm->d_ins_rowIndices, ins_rowIndices, sizeof(int)*(1+(int)ceild(numRows,(BLOCKSIZE/HALFWARP))), cudaMemcpyHostToDevice);
    cudaMemcpy(spm->d_ins_inputList, ins_inputList, sizeof(int)*ins_inputListCount, cudaMemcpyHostToDevice);
  #endif
}

