/*
 * IBM Sparse Matrix-Vector Multiplication Toolkit for Graphics Processing Units
 * (c) Copyright IBM Corp. 2008, 2009.  All Rights Reserved.
 */ 

#ifndef _SPMV_KERNEL_H_
#define _SPMV_KERNEL_H_

#if CACHE
    texture<float,1> tex_y_float;
#endif


////////////////////////////////////////////////////////////////////////////////
// SpMV Kernel Device Code
////////////////////////////////////////////////////////////////////////////////
__global__ void
SpMV(float *x, const float *val, const unsigned int *rowIndices,
               const unsigned int *indices, const float *y, const unsigned int numRows, 
	       const unsigned int numCols, const unsigned int numNonZeroElements)
{
    unsigned int tid = threadIdx.y;
    unsigned int bid = blockIdx.y;
    float t=0;
#if C_GLOBAL_OPT
    unsigned int ind2Dx = tid%HALFWARP;
    unsigned int ind2Dy = tid/HALFWARP;
    unsigned int ub,lb;
    unsigned int myblock = bid * (BLOCKSIZE/HALFWARP);
    unsigned int myi = myblock + ind2Dy;
    __shared__ int rowInd[(BLOCKSIZE/HALFWARP)+1];
    __shared__ float tempProd[(BLOCKSIZE/HALFWARP)][HALFWARP+PAD];

    if ((tid <= ((BLOCKSIZE/HALFWARP))) && (myi < numRows))
        rowInd[tid] = rowIndices[myblock+tid];
    __syncthreads();
    t=0;
    unsigned int rowStartNZ = rowInd[ind2Dy];
    unsigned int rowStartRes = rowStartNZ%HALFWARP;
    unsigned int rowStartAlignNZ = rowStartNZ + HALFWARP - rowStartRes;
    lb = rowStartAlignNZ+ind2Dx;
    ub = rowInd[ind2Dy+1];
    unsigned int j ;
    if (myi < numRows) {
        j = rowStartAlignNZ - HALFWARP + ind2Dx;
        if ( (j >= rowStartNZ) && (j<ub) ) { 
            unsigned int ind = indices[j];
	#if CACHE
	    float yval = tex1Dfetch(tex_y_float, ind);
	#else
	    float yval = y[ind];
	#endif
            t += val[j] * yval;
	}
        for (j=lb; j<ub; j+=HALFWARP) {
            unsigned int ind = indices[j];
	#if CACHE
	    float yval = tex1Dfetch(tex_y_float, ind);
	#else
	    float yval = y[ind];
	#endif
            t += val[j] * yval;
        }
        tempProd[ind2Dy][ind2Dx] = t;
    }
    __syncthreads();
#if 0
    if ((ind2Dx == 0) && (myi < numRows)) {
        t=0;
        for (unsigned int k = 0; k<HALFWARP; k++) {
            t += tempProd[ind2Dy][k];
        }
        x[myi] = t;
    }
    __syncthreads();
#endif
#if 1
    // Works for HALFWARP=16
    if ((ind2Dx == 0) && (myi < numRows)) {
        t = tempProd[ind2Dy][0] + tempProd[ind2Dy][1] + tempProd[ind2Dy][2] + tempProd[ind2Dy][3] +\
          tempProd[ind2Dy][4] + tempProd[ind2Dy][5] + tempProd[ind2Dy][6] + tempProd[ind2Dy][7] +\
          tempProd[ind2Dy][8] + tempProd[ind2Dy][9] + tempProd[ind2Dy][10]+ tempProd[ind2Dy][11]+\
          tempProd[ind2Dy][12]+ tempProd[ind2Dy][13]+ tempProd[ind2Dy][14]+ tempProd[ind2Dy][15];
        x[myi] = t;
    }
    __syncthreads();
#endif
#endif
#if NEW_GLOBAL_OPT
    unsigned int ind2Dx = tid%HALFWARP;
    unsigned int ind2Dy = tid/HALFWARP;
    unsigned int ub,lb;
    unsigned int blockStart = bid * (BLOCKSIZE/HALFWARP);
    unsigned int blockEnd = min(blockStart+(BLOCKSIZE/HALFWARP)-1,numRows-1);
    unsigned int myi = blockStart + ind2Dy;

    __shared__ int rowInd[(BLOCKSIZE/HALFWARP)+1];
    __shared__ float tempProd[(BLOCKSIZE/HALFWARP)][(BLOCKSIZE/HALFWARP)][HALFWARP+PAD];

    if ((tid <= ((BLOCKSIZE/HALFWARP))) && ((blockStart+tid) <= numRows))
        rowInd[tid] = rowIndices[blockStart+tid];
    __syncthreads();

    t=0;
    lb = rowInd[0]+tid;
    ub = rowInd[blockEnd-blockStart+1];
    unsigned int curr_i=0;
    for (int p=0;p<(BLOCKSIZE/HALFWARP);p++)
        tempProd[p][ind2Dy][ind2Dx] = 0;
    __syncthreads();
    for (unsigned int j=lb; j<ub; j+=NUMTHREADS) {
    	for (int p=curr_i;p<(blockEnd-blockStart+1);p++) {
	    if ( (j >= rowInd[p]) && (j < rowInd[p+1]) ) {
	        curr_i = p;
	        break;
	    }
	}
        unsigned int ind = indices[j];
	#if CACHE
	    float yval = tex1Dfetch(tex_y_float, ind);
	#else
	    float yval = y[ind];
	#endif
        t = val[j] * yval;
        tempProd[curr_i][ind2Dy][ind2Dx] += t;
    }
    __syncthreads();
#if 0
    if ((ind2Dx == 0) && (myi < numRows)) {
        t=0;
        for (unsigned int k = 0; k<(BLOCKSIZE/HALFWARP); k++) {
	    for (unsigned int l = 0; l<HALFWARP; l++) 
		t += tempProd[ind2Dy][k][l];
        }
        x[myi] = t;
    }
#endif
#if 1
    if ((ind2Dx == 0) && (myi < numRows)) {
	t=0;
	for (unsigned int k = 0; k<(BLOCKSIZE/HALFWARP); k++) {
            t += tempProd[ind2Dy][k][0] + tempProd[ind2Dy][k][1] + tempProd[ind2Dy][k][2] + tempProd[ind2Dy][k][3] +\
              tempProd[ind2Dy][k][4] + tempProd[ind2Dy][k][5] + tempProd[ind2Dy][k][6] + tempProd[ind2Dy][k][7] +\
              tempProd[ind2Dy][k][8] + tempProd[ind2Dy][k][9] + tempProd[ind2Dy][k][10]+ tempProd[ind2Dy][k][11]+\
              tempProd[ind2Dy][k][12]+ tempProd[ind2Dy][k][13]+ tempProd[ind2Dy][k][14]+ tempProd[ind2Dy][k][15];
	}
        x[myi] = t;
    }
#endif
     __syncthreads();
#endif
#if GLOBAL_OPT
    unsigned int ind2Dx = tid%HALFWARP;
    unsigned int ind2Dy = tid/HALFWARP;
    unsigned int ub,lb;
    unsigned int myblock = bid * (BLOCKSIZE/HALFWARP);
    unsigned int myi = myblock + ind2Dy;
    __shared__ int rowInd[(BLOCKSIZE/HALFWARP)+1];
    __shared__ float tempProd[(BLOCKSIZE/HALFWARP)][HALFWARP+PAD];
	if ((tid <= ((BLOCKSIZE/HALFWARP))) && (myi < numRows)) 
            rowInd[tid] = rowIndices[myblock+tid];
	__syncthreads();
	t=0;
	lb = rowInd[ind2Dy]+ind2Dx;
	ub = rowInd[ind2Dy+1];
	if (myi < numRows) {
            for (unsigned int j=lb; j<ub; j+=HALFWARP) {
       	    	unsigned int ind = indices[j];
	     #if CACHE
	    	float yval = tex1Dfetch(tex_y_float, ind);
	     #else
	    	float yval = y[ind];
	     #endif
		t += val[j] * yval;
	    }
	    tempProd[ind2Dy][ind2Dx] = t;
	}
	//__syncthreads();
#if 0 
	if ((ind2Dx == 0) && (myi < numRows)) {
	    t=0;
	    for (unsigned int k = 0; k<HALFWARP; k++) {
		t += tempProd[ind2Dy][k];
	    }
       	    x[myi] = t;
	}
	__syncthreads();
#endif
#if 0
	// Works for HALFWARP=8
	if ((ind2Dx == 0) && (myi < numRows)) {
            t = tempProd[ind2Dy][0] + tempProd[ind2Dy][1] + tempProd[ind2Dy][2] + tempProd[ind2Dy][3] +\
              tempProd[ind2Dy][4] + tempProd[ind2Dy][5] + tempProd[ind2Dy][6] + tempProd[ind2Dy][7];
       	    x[myi] = t;
	}
#endif
#if 1
	// Works for HALFWARP=16
	if ((ind2Dx == 0) && (myi < numRows)) {
            t = tempProd[ind2Dy][0] + tempProd[ind2Dy][1] + tempProd[ind2Dy][2] + tempProd[ind2Dy][3] +\
              tempProd[ind2Dy][4] + tempProd[ind2Dy][5] + tempProd[ind2Dy][6] + tempProd[ind2Dy][7] +\
              tempProd[ind2Dy][8] + tempProd[ind2Dy][9] + tempProd[ind2Dy][10]+ tempProd[ind2Dy][11]+\
              tempProd[ind2Dy][12]+ tempProd[ind2Dy][13]+ tempProd[ind2Dy][14]+ tempProd[ind2Dy][15];
       	    x[myi] = t;
	}
	//__syncthreads();
#endif
#if 0
	// Works for HALFWARP=16/32
	if (myi < numRows) {
	    //if (ind2Dx < 16) tempProd[ind2Dy][ind2Dx] += tempProd[ind2Dy][ind2Dx+16];
	    if (ind2Dx < 8) tempProd[ind2Dy][ind2Dx] += tempProd[ind2Dy][ind2Dx+8];
	    if (ind2Dx < 4) tempProd[ind2Dy][ind2Dx] += tempProd[ind2Dy][ind2Dx+4];
	    if (ind2Dx < 2) tempProd[ind2Dy][ind2Dx] += tempProd[ind2Dy][ind2Dx+2];
	    if (ind2Dx < 1) x[myi]= tempProd[ind2Dy][ind2Dx] +tempProd[ind2Dy][ind2Dx+1];
	}
	__syncthreads();
#endif
#if 0
        // Works for HALFWARP=16 & 32
        if (!(ind2Dx % 4) && (myi < numRows)) {
            tempProd[ind2Dy][ind2Dx] += tempProd[ind2Dy][ind2Dx+1] + tempProd[ind2Dy][ind2Dx+2] + tempProd[ind2Dy][ind2Dx+3];
        }
        __syncthreads();
        if ((ind2Dx == 0) && (myi < numRows)) {
        #if 1 // for halfwarp 16
            t = tempProd[ind2Dy][0] + tempProd[ind2Dy][4] + tempProd[ind2Dy][8] + tempProd[ind2Dy][12];
        #else // for halfwarp 32
            t = tempProd[ind2Dy][0] + tempProd[ind2Dy][4] + tempProd[ind2Dy][8] + tempProd[ind2Dy][12]+\
               +tempProd[ind2Dy][16] + tempProd[ind2Dy][20] + tempProd[ind2Dy][24] + tempProd[ind2Dy][28];
        #endif
            x[myi] = t;
        }
        __syncthreads();
#endif
#endif
#if GLOBAL
    unsigned int myi = bid * BLOCKSIZE + tid;
    if (myi < numRows) {
    	unsigned int lb = rowIndices[myi];
    	unsigned int ub = rowIndices[myi+1];
    	for (unsigned int j=lb; j<ub; j++) {
	    unsigned int ind = indices[j];
        #if CACHE
            float yval = tex1Dfetch(tex_y_float, ind);
        #else
            float yval = y[ind];
        #endif
            t += val[j] * yval;
    	}
    	x[myi] = t;
    }
    __syncthreads();
#endif
#if SHARED_RI
    unsigned int myi = bid * BLOCKSIZE + tid;
    __shared__ int rowInd[NUMTHREADS+1];
    if (myi < numRows) {
    	rowInd[tid] = rowIndices[myi];
    }
    __syncthreads();
    if ( (myi < numRows) && ((tid == NUMTHREADS-1) || (myi == numRows-1) ) ) {
	rowInd[tid+1] = rowIndices[myi+1];
    }
    __syncthreads();
    if (myi < numRows) {
    	unsigned int lb = rowInd[tid];
    	unsigned int ub = rowInd[tid+1];
    	for (unsigned int j=lb; j<ub; j++) {
            unsigned int ind = indices[j];
        #if CACHE
            float yval = tex1Dfetch(tex_y_float, ind);
        #else
            float yval = y[ind];
        #endif
            t += val[j] * yval;
    	}
    	x[myi] = t;
    }
    __syncthreads();
#endif
#if GLOBAL_SHARED_OPT
    unsigned int ind2Dx = tid%HALFWARP;
    unsigned int ind2Dy = tid/HALFWARP;
    unsigned int myblock = bid * (BLOCKSIZE/HALFWARP);
    unsigned int myi = myblock + ind2Dy;
    __shared__ int rowInd[(BLOCKSIZE/HALFWARP)+1];
    __shared__ float tempProd[(BLOCKSIZE/HALFWARP)][HALFWARP+PAD];
    __shared__ float ys[NUMTHREADS];
        if ((tid <= ((BLOCKSIZE/HALFWARP))) && (myi < numRows))
            rowInd[tid] = rowIndices[myblock+tid];
        __syncthreads();
        t=0;
        unsigned int ind, lb, ub, j;
	if (myi < numRows) {
	    lb = rowInd[ind2Dy]+ind2Dx;
	    ub = rowInd[ind2Dy+1];
	    j = lb;
	    if (j<ub)	ind = indices[j];
	    else ind = numCols+BLOCKSIZE;
	}
        __syncthreads();
    	for (unsigned int k=0; k<numCols; k+=BLOCKSIZE) {
            __syncthreads();
            if ( (k+tid) < numCols)
                ys[tid] = y[k+tid];
            __syncthreads();

            if (myi < numRows) {
	    #if 0		
            	while ( ((j+HALFWARP) < ub) && ( ind < (k+BLOCKSIZE) )  ) {
                    t += val[j] * ys[ind-k];
                    j+=HALFWARP;
                    ind = indices[j];
            	}
	    #endif
	    #if 1
            	while ( ind < (k+BLOCKSIZE) ) {
                    t += val[j] * ys[ind-k];
                    j+=HALFWARP;
		    if (j < ub) ind = indices[j];
		    else { ind = numCols+BLOCKSIZE; break; }
            	}
	    #endif
            }
	}
     #if 0
	if ( (myi < numRows) && (j<ub) )
            tempProd[ind2Dy][ind2Dx] = t + val[j] * y[ind];
	else
	    tempProd[ind2Dy][ind2Dx] = t;
     #endif
     #if 1
        tempProd[ind2Dy][ind2Dx] = t;
     #endif
	__syncthreads();
#if 0
        if ((ind2Dx == 0) && (myi < numRows)) {
            t=0;
            for (unsigned int k = 0; k<HALFWARP; k++) {
                t += tempProd[ind2Dy][k];
            }
            x[myi] = t;
        }
        __syncthreads();
#endif
#if 0
        // Works for HALFWARP=16
        if ((ind2Dx == 0) && (myi < numRows)) {
            t = tempProd[ind2Dy][0] + tempProd[ind2Dy][1] + tempProd[ind2Dy][2] + tempProd[ind2Dy][3] +\
              tempProd[ind2Dy][4] + tempProd[ind2Dy][5] + tempProd[ind2Dy][6] + tempProd[ind2Dy][7] +\
              tempProd[ind2Dy][8] + tempProd[ind2Dy][9] + tempProd[ind2Dy][10]+ tempProd[ind2Dy][11]+\
              tempProd[ind2Dy][12]+ tempProd[ind2Dy][13]+ tempProd[ind2Dy][14]+ tempProd[ind2Dy][15];
            x[myi] = t;
        }
        __syncthreads();
#endif
#if 1
        // Works for HALFWARP=16 & 32
        if (!(ind2Dx % 4) && (myi < numRows)) {
            tempProd[ind2Dy][ind2Dx] += tempProd[ind2Dy][ind2Dx+1] + tempProd[ind2Dy][ind2Dx+2] + tempProd[ind2Dy][ind2Dx+3];
        }
        __syncthreads();
        if ((ind2Dx == 0) && (myi < numRows)) {
        #if 1 // for halfwarp 16
            t = tempProd[ind2Dy][0] + tempProd[ind2Dy][4] + tempProd[ind2Dy][8] + tempProd[ind2Dy][12];
        #else // for halfwarp 32
            t = tempProd[ind2Dy][0] + tempProd[ind2Dy][4] + tempProd[ind2Dy][8] + tempProd[ind2Dy][12]+\
               +tempProd[ind2Dy][16] + tempProd[ind2Dy][20] + tempProd[ind2Dy][24] + tempProd[ind2Dy][28];
        #endif
            x[myi] = t;
        }
        __syncthreads();
#endif
#endif

}

__global__ void
SpMV_withInspectInput(float *x, const float *val, const unsigned int *rowIndices,
                const unsigned int *indices, const float *y, const unsigned int numRows,
                const unsigned int numCols, const unsigned int numNonZeroElements,
                const unsigned int *ins_rowIndices, const unsigned int *ins_indices, const unsigned int *ins_inputList)
{
#if CACHE
    unsigned int tid = threadIdx.y;
    unsigned int bid = blockIdx.y;
    float t=0;
    unsigned int ind2Dx = tid%HALFWARP;
    unsigned int ind2Dy = tid/HALFWARP;
    unsigned int myblock = bid * (BLOCKSIZE/HALFWARP);
    unsigned int myi = myblock + ind2Dy;
    __shared__ int rowInd[(BLOCKSIZE/HALFWARP)+1];
    __shared__ float tempProd[(BLOCKSIZE/HALFWARP)][HALFWARP+PAD];
    __shared__ float ys[INSPECT_INPUT_MAX];
    __shared__ int ins_rowInd[2];
    __shared__ int ins_Ind[BLOCKSIZE]; // Have to fix
    __shared__ int ins_inpStartVal[BLOCKSIZE]; // Have to fix

    if ((tid <= ((BLOCKSIZE/HALFWARP))) && (myi < numRows))
        rowInd[tid] = rowIndices[myblock+tid];
    __syncthreads();
    if ((ind2Dx < 2) && (ind2Dy == 0) && (myi < numRows))
        ins_rowInd[ind2Dx]=ins_rowIndices[bid+ind2Dx];
    __syncthreads();

    t=0;
    unsigned int lb = rowInd[ind2Dy]+ind2Dx;
    unsigned int ub = rowInd[ind2Dy+1];
    unsigned int ktlb = ins_rowInd[0];
    unsigned int ktub = ins_rowInd[1];
    if (tid <= (ktub-ktlb)) ins_Ind[tid]=ins_indices[tid+ktlb];
    __syncthreads();
    if (tid < (ktub-ktlb)) {
	unsigned int is = ins_Ind[tid];
	ins_inpStartVal[tid]=ins_inputList[is];
    }
    __syncthreads();
    unsigned int kt=ktlb;
    unsigned int j=lb;
    for (;kt<ktub;kt++) {
	unsigned int startVal = ins_inpStartVal[kt-ktlb];
        if (startVal != numCols) {
                
            unsigned int is = ins_Ind[kt-ktlb];
            unsigned int ie = ins_Ind[kt-ktlb+1];
            for (unsigned int iL=is+tid;iL<ie;iL+=NUMTHREADS) {
		unsigned int currInd = ins_inputList[iL];
		ys[currInd-startVal] = tex1Dfetch(tex_y_float, currInd);	
	    }
	    __syncthreads();
    	    if (myi < numRows && j<ub) {
		unsigned int ind = indices[j];
		t += val[j] * ys[ind-startVal];
		j+=HALFWARP;
	    }
  	}
	else {
    	    if (myi < numRows && j<ub) {
                unsigned int ind = indices[j];
                float yval = tex1Dfetch(tex_y_float, ind);
                t += val[j] * yval;
		j+=HALFWARP;
	    }
	}
    }
    tempProd[ind2Dy][ind2Dx] = t;
   // __syncthreads();
  #if 0
    if ((ind2Dx == 0) && (myi < numRows)) {
        t=0;
        for (unsigned int k = 0; k<HALFWARP; k++) {
            t += tempProd[ind2Dy][k];
        }
        x[myi] = t;
    }
    __syncthreads();
  #endif
  #if 1
    // Works for HALFWARP=16
    if ((ind2Dx == 0) && (myi < numRows)) {
        t = tempProd[ind2Dy][0] + tempProd[ind2Dy][1] + tempProd[ind2Dy][2] + tempProd[ind2Dy][3] +\
          tempProd[ind2Dy][4] + tempProd[ind2Dy][5] + tempProd[ind2Dy][6] + tempProd[ind2Dy][7] +\
          tempProd[ind2Dy][8] + tempProd[ind2Dy][9] + tempProd[ind2Dy][10]+ tempProd[ind2Dy][11]+\
          tempProd[ind2Dy][12]+ tempProd[ind2Dy][13]+ tempProd[ind2Dy][14]+ tempProd[ind2Dy][15];
        x[myi] = t;
    }
   // __syncthreads();
  #endif
#endif
}


__global__ void
SpMV_withInspect(float *x, const float *val, const unsigned int *rowIndices,
                const unsigned int *indices, const float *y, const unsigned int numRows,
                const unsigned int numCols, const unsigned int numNonZeroElements,
                const unsigned int *ins_rowIndices, const unsigned int *ins_indices)
{
#if C_GLOBAL_OPT
#if 0 
    unsigned int tid = threadIdx.y;
    unsigned int bid = blockIdx.y;
    float t=0;
    unsigned int ind2Dx = tid%HALFWARP;
    unsigned int ind2Dy = tid/HALFWARP;
    unsigned int myblock = bid * (BLOCKSIZE/HALFWARP);
    unsigned int myi = myblock + ind2Dy;
    __shared__ int rowInd[(BLOCKSIZE/HALFWARP)+1];
    __shared__ float tempProd[(BLOCKSIZE/HALFWARP)][HALFWARP+PAD];
    __shared__ float ys[INSPECT_BLOCK_c];
    __shared__ int ins_rowInd[2];
    __shared__ int ins_Ind;

    if ((tid <= ((BLOCKSIZE/HALFWARP))) && (myi < numRows))
        rowInd[tid] = rowIndices[myblock+tid];
    __syncthreads();
    if ((ind2Dx < 2) && (ind2Dy == 0) && (myi < numRows))
        ins_rowInd[ind2Dx]=ins_rowIndices[bid+ind2Dx];
    __syncthreads();

    t=0;
    unsigned int ind, lb, ub, j;
    float valS;
    if (myi < numRows) {
        lb = rowInd[ind2Dy]+ind2Dx;
        ub = rowInd[ind2Dy+1];
        j = lb;
        ind = indices[j];
	valS = val[j]; 
    }
    __syncthreads();
    unsigned int ktlb = ins_rowInd[0];
    unsigned int ktub = ins_rowInd[1];
    for (unsigned int kt=ktlb; kt<ktub; kt++) {
        __syncthreads();
        if (tid==0) ins_Ind=ins_indices[kt];
        __syncthreads();
    #if VAR_BLOCK
        unsigned int k = ins_Ind; // In case of var_block, ins_indices 'll have original column index
    #else
        unsigned int k = ins_Ind*INSPECT_BLOCK_c;
    #endif
        if ( tid < min(INSPECT_BLOCK_c, numCols-k) )
	  #if CACHE
            ys[tid] = tex1Dfetch(tex_y_float, k+tid);
	  #else
            ys[tid] = y[k+tid];
	  #endif
        __syncthreads();
        if (myi < numRows) {
            while (ind < (k+INSPECT_BLOCK_c)) {
                t += valS * ys[ind-k];
                j+=HALFWARP;
                if (j < ub) { ind = indices[j]; valS = val[j]; }
                else { ind = 2*numCols;  }
            }
        }
    }
    tempProd[ind2Dy][ind2Dx] = t;
    __syncthreads();
  #if 0
    if ((ind2Dx == 0) && (myi < numRows)) {
        t=0;
        for (unsigned int k = 0; k<HALFWARP; k++) {
            t += tempProd[ind2Dy][k];
        }
        x[myi] = t;
    }
    __syncthreads();
  #endif
  #if 1
    // Works for HALFWARP=16
    if ((ind2Dx == 0) && (myi < numRows)) {
        t = tempProd[ind2Dy][0] + tempProd[ind2Dy][1] + tempProd[ind2Dy][2] + tempProd[ind2Dy][3] +\
          tempProd[ind2Dy][4] + tempProd[ind2Dy][5] + tempProd[ind2Dy][6] + tempProd[ind2Dy][7] +\
          tempProd[ind2Dy][8] + tempProd[ind2Dy][9] + tempProd[ind2Dy][10]+ tempProd[ind2Dy][11]+\
          tempProd[ind2Dy][12]+ tempProd[ind2Dy][13]+ tempProd[ind2Dy][14]+ tempProd[ind2Dy][15];
        x[myi] = t;
    }
    __syncthreads();
  #endif
#endif
#if 1 
    unsigned int tid = threadIdx.y;
    unsigned int bid = blockIdx.y;
    float t=0;
    unsigned int ind2Dx = tid%HALFWARP;
    unsigned int ind2Dy = tid/HALFWARP;
    unsigned int myblock = bid * (BLOCKSIZE/HALFWARP);
    unsigned int myi = myblock + ind2Dy;
    __shared__ int rowInd[(BLOCKSIZE/HALFWARP)+1];
    __shared__ float tempProd[(BLOCKSIZE/HALFWARP)][HALFWARP+PAD];
    __shared__ float ys[INSPECT_BLOCK_c];
    __shared__ int ins_rowInd[2];
    //__shared__ int ins_Ind;
    __shared__ int ins_Ind[BLOCKSIZE];

    if ((tid <= ((BLOCKSIZE/HALFWARP))) && (myi < numRows))
        rowInd[tid] = rowIndices[myblock+tid];
    //__syncthreads();
    if ((ind2Dx < 2) && (ind2Dy == 0) && (myi < numRows))
        ins_rowInd[ind2Dx]=ins_rowIndices[bid+ind2Dx];
    __syncthreads();

    t=0;
    unsigned int ind, lb, ub, j;
    float valS;
    if (myi < numRows) {
        lb = rowInd[ind2Dy]+ind2Dx;
        ub = rowInd[ind2Dy+1];
        j = lb;
        ind = indices[j];
	valS = val[j]; 
    }
    unsigned int ktlb = ins_rowInd[0];
    unsigned int ktub = ins_rowInd[1];
    if (tid < (ktub-ktlb)) ins_Ind[tid]=ins_indices[tid+ktlb];
    __syncthreads();

    for (unsigned int kt=ktlb; kt<ktub; kt++) {
      #if VAR_BLOCK
        unsigned int k = ins_Ind[kt-ktlb]; // In case of var_block, ins_indices 'll have original column index
      #else
        unsigned int k = ins_Ind[kt-ktlb]*INSPECT_BLOCK_c;
      #endif
        if ( tid < min(INSPECT_BLOCK_c, numCols-k) )
	  #if CACHE
            ys[tid] = tex1Dfetch(tex_y_float, k+tid);
	  #else
            ys[tid] = y[k+tid];
	  #endif
        __syncthreads();
        if (myi < numRows) {
            while (ind < (k+INSPECT_BLOCK_c)) {
                t += valS * ys[ind-k];
                j+=HALFWARP;
                if (j < ub) { ind = indices[j]; valS = val[j]; }
                else { ind = 2*numCols;  }
            }
        }
        //__syncthreads();
    }

    tempProd[ind2Dy][ind2Dx] = t;
    //__syncthreads();
  #if 0
    if ((ind2Dx == 0) && (myi < numRows)) {
        t=0;
        for (unsigned int k = 0; k<HALFWARP; k++) {
            t += tempProd[ind2Dy][k];
        }
        x[myi] = t;
    }
    __syncthreads();
  #endif
  #if 1
    // Works for HALFWARP=16
    if ((ind2Dx == 0) && (myi < numRows)) {
        t = tempProd[ind2Dy][0] + tempProd[ind2Dy][1] + tempProd[ind2Dy][2] + tempProd[ind2Dy][3] +\
          tempProd[ind2Dy][4] + tempProd[ind2Dy][5] + tempProd[ind2Dy][6] + tempProd[ind2Dy][7] +\
          tempProd[ind2Dy][8] + tempProd[ind2Dy][9] + tempProd[ind2Dy][10]+ tempProd[ind2Dy][11]+\
          tempProd[ind2Dy][12]+ tempProd[ind2Dy][13]+ tempProd[ind2Dy][14]+ tempProd[ind2Dy][15];
        x[myi] = t;
    }
    //__syncthreads();
  #endif
  #if 0
    // Works for HALFWARP=8
    if ((ind2Dx == 0) && (myi < numRows)) {
        t = tempProd[ind2Dy][0] + tempProd[ind2Dy][1] + tempProd[ind2Dy][2] + tempProd[ind2Dy][3] +\
          tempProd[ind2Dy][4] + tempProd[ind2Dy][5] + tempProd[ind2Dy][6] + tempProd[ind2Dy][7];
        x[myi] = t;
    }
    //__syncthreads();
  #endif
#endif
#if 0 
    unsigned int tid = threadIdx.y;
    unsigned int bid = blockIdx.y;
    float t=0;
    unsigned int ind2Dx = tid%HALFWARP;
    unsigned int ind2Dy = tid/HALFWARP;
    unsigned int myblock = bid * (BLOCKSIZE/HALFWARP);
    unsigned int myi = myblock + ind2Dy;
    __shared__ int rowInd[(BLOCKSIZE/HALFWARP)+1];
    __shared__ float tempProd[(BLOCKSIZE/HALFWARP)][HALFWARP+PAD];
    __shared__ float ys[INSPECT_BLOCK_c];
    __shared__ int ins_rowInd[2];
    //__shared__ int ins_Ind;
    __shared__ int ins_Ind[BLOCKSIZE];

    if ((tid <= ((BLOCKSIZE/HALFWARP))) && (myi < numRows))
        rowInd[tid] = rowIndices[myblock+tid];
    __syncthreads();
    if ((ind2Dx < 2) && (ind2Dy == 0) && (myi < numRows))
        ins_rowInd[ind2Dx]=ins_rowIndices[bid+ind2Dx];
    __syncthreads();

    t=0;
    unsigned int ind, lb, ub, j;
    float valS;
    if (myi < numRows) {
        lb = rowInd[ind2Dy]+ind2Dx;
        ub = rowInd[ind2Dy+1];
        j = lb;
        ind = indices[j];
	valS = val[j]; 
    }
    unsigned int ktlb = ins_rowInd[0];
    unsigned int ktub = ins_rowInd[1];
    if (tid < (ktub-ktlb)) ins_Ind[tid]=ins_indices[tid+ktlb];
    __syncthreads();

    unsigned int kt=ktlb;
    if((ktub-ktlb)%2) {
      #if VAR_BLOCK
        unsigned int k1 = ins_Ind[kt-ktlb]; // In case of var_block, ins_indices 'll have original column index
      #else
        unsigned int k1 = ins_Ind[kt-ktlb]*INSPECT_BLOCK_c;
      #endif
        if ( tid < min(INSPECT_BLOCK_c, numCols-k1) )
            ys[tid] = y[k1+tid];
        __syncthreads();
        if (myi < numRows) {
            while (ind < (k1+INSPECT_BLOCK_c)) {
                t += valS * ys[ind-k1];
                j+=HALFWARP;
                if (j < ub) { ind = indices[j]; valS = val[j]; }
                else { ind = 2*numCols;  }
            }
        }
        //__syncthreads();
	kt++;
    }
    for (; kt<ktub; kt+=2) {
      #if VAR_BLOCK
        unsigned int k1 = ins_Ind[kt-ktlb]; // In case of var_block, ins_indices 'll have original column index
        unsigned int k2 = ins_Ind[kt-ktlb+1]; // In case of var_block, ins_indices 'll have original column index
      #else
        unsigned int k1 = ins_Ind[kt-ktlb]*INSPECT_BLOCK_c;
        unsigned int k2 = ins_Ind[kt-ktlb+1]*INSPECT_BLOCK_c;
      #endif
        if ( tid < min(INSPECT_BLOCK_c, numCols-k1) )
            ys[tid] = y[k1+tid];
        __syncthreads();
        if (myi < numRows) {
            while (ind < (k1+INSPECT_BLOCK_c)) {
                t += valS * ys[ind-k1];
                j+=HALFWARP;
                if (j < ub) { ind = indices[j]; valS = val[j]; }
                else { ind = 2*numCols;  }
            }
        }
        //__syncthreads();
        if ( tid < min(INSPECT_BLOCK_c, numCols-k2) )
            ys[tid] = y[k2+tid];
        __syncthreads();
        if (myi < numRows) {
            while (ind < (k2+INSPECT_BLOCK_c)) {
                t += valS * ys[ind-k2];
                j+=HALFWARP;
                if (j < ub) { ind = indices[j]; valS = val[j]; }
                else { ind = 2*numCols;  }
            }
        }
        //__syncthreads();
    }

    tempProd[ind2Dy][ind2Dx] = t;
    //__syncthreads();
  #if 0
    if ((ind2Dx == 0) && (myi < numRows)) {
        t=0;
        for (unsigned int k = 0; k<HALFWARP; k++) {
            t += tempProd[ind2Dy][k];
        }
        x[myi] = t;
    }
    __syncthreads();
  #endif
  #if 1
    // Works for HALFWARP=16
    if ((ind2Dx == 0) && (myi < numRows)) {
        t = tempProd[ind2Dy][0] + tempProd[ind2Dy][1] + tempProd[ind2Dy][2] + tempProd[ind2Dy][3] +\
          tempProd[ind2Dy][4] + tempProd[ind2Dy][5] + tempProd[ind2Dy][6] + tempProd[ind2Dy][7] +\
          tempProd[ind2Dy][8] + tempProd[ind2Dy][9] + tempProd[ind2Dy][10]+ tempProd[ind2Dy][11]+\
          tempProd[ind2Dy][12]+ tempProd[ind2Dy][13]+ tempProd[ind2Dy][14]+ tempProd[ind2Dy][15];
        x[myi] = t;
    }
    //__syncthreads();
  #endif
  #if 0
    // Works for HALFWARP=8
    if ((ind2Dx == 0) && (myi < numRows)) {
        t = tempProd[ind2Dy][0] + tempProd[ind2Dy][1] + tempProd[ind2Dy][2] + tempProd[ind2Dy][3] +\
          tempProd[ind2Dy][4] + tempProd[ind2Dy][5] + tempProd[ind2Dy][6] + tempProd[ind2Dy][7];
        x[myi] = t;
    }
    //__syncthreads();
  #endif
#endif
#if 0
    unsigned int tid = threadIdx.y;
    unsigned int bid = blockIdx.y;
    float t=0;
    unsigned int ind2Dx = tid%HALFWARP;
    unsigned int ind2Dy = tid/HALFWARP;
    unsigned int myblock = bid * (BLOCKSIZE/HALFWARP);
    unsigned int myi = myblock + ind2Dy;
    __shared__ int rowInd[(BLOCKSIZE/HALFWARP)+1];
    __shared__ float tempProd[(BLOCKSIZE/HALFWARP)][HALFWARP+PAD];
    __shared__ float ys[INSPECT_BLOCK_c];
    __shared__ int ins_rowInd[2];
    __shared__ int ins_Ind;

    if ((tid <= ((BLOCKSIZE/HALFWARP))) && (myi < numRows))
        rowInd[tid] = rowIndices[myblock+tid];
    //__syncthreads();
    if ((ind2Dx < 2) && (ind2Dy == 0) && (myi < numRows))
        ins_rowInd[ind2Dx]=ins_rowIndices[bid+ind2Dx];
    __syncthreads();

    t=0;
    unsigned int ind, lb, ub, j;
    float valS;
    if (myi < numRows) {
        lb = rowInd[ind2Dy]+ind2Dx;
        ub = rowInd[ind2Dy+1];
        j = lb;
        //ind = indices[j];
        //valS = val[j];
    }
    __syncthreads();
    unsigned int ktlb = ins_rowInd[0];
    unsigned int ktub = ins_rowInd[1];
    for (unsigned int kt=ktlb; kt<ktub; kt++) {
        //__syncthreads();
        if (tid==0) ins_Ind=ins_indices[kt];
        __syncthreads();
    #if VAR_BLOCK
        unsigned int k = ins_Ind; // In case of var_block, ins_indices 'll have original column index
    #else
        unsigned int k = ins_Ind*INSPECT_BLOCK_c;
    #endif
        if ( tid < min(INSPECT_BLOCK_c, numCols-k) )
            ys[tid] = y[k+tid];
        __syncthreads();
        if (myi < numRows) {
            //while (ind < (k+INSPECT_BLOCK_c)) {
                ind = indices[j]; valS = val[j];
		t += valS * ys[ind-k];
                j+=HALFWARP;
                //if (j < ub) { ind = indices[j]; valS = val[j]; }
                //else { ind = 2*numCols;  }
            //}
        }
    }
    tempProd[ind2Dy][ind2Dx] = t;
    //__syncthreads();
  #if 0
    if ((ind2Dx == 0) && (myi < numRows)) {
        t=0;
        for (unsigned int k = 0; k<HALFWARP; k++) {
            t += tempProd[ind2Dy][k];
        }
        x[myi] = t;
    }
    __syncthreads();
  #endif
  #if 1
    // Works for HALFWARP=16
    if ((ind2Dx == 0) && (myi < numRows)) {
        t = tempProd[ind2Dy][0] + tempProd[ind2Dy][1] + tempProd[ind2Dy][2] + tempProd[ind2Dy][3] +\
          tempProd[ind2Dy][4] + tempProd[ind2Dy][5] + tempProd[ind2Dy][6] + tempProd[ind2Dy][7] +\
          tempProd[ind2Dy][8] + tempProd[ind2Dy][9] + tempProd[ind2Dy][10]+ tempProd[ind2Dy][11]+\
          tempProd[ind2Dy][12]+ tempProd[ind2Dy][13]+ tempProd[ind2Dy][14]+ tempProd[ind2Dy][15];
        x[myi] = t;
    }
  #endif
  #if 0
    // Works for HALFWARP=8
    if ((ind2Dx == 0) && (myi < numRows)) {
        t = tempProd[ind2Dy][0] + tempProd[ind2Dy][1] + tempProd[ind2Dy][2] + tempProd[ind2Dy][3] +\
          tempProd[ind2Dy][4] + tempProd[ind2Dy][5] + tempProd[ind2Dy][6] + tempProd[ind2Dy][7];
        x[myi] = t;
    }
    //__syncthreads();
  #endif
#endif
#else
    unsigned int tid = threadIdx.y;
    unsigned int bid = blockIdx.y;
    float t=0;
    unsigned int ind2Dx = tid%HALFWARP;
    unsigned int ind2Dy = tid/HALFWARP;
    unsigned int myblock = bid * (BLOCKSIZE/HALFWARP);
    unsigned int myi = myblock + ind2Dy;
    __shared__ int rowInd[(BLOCKSIZE/HALFWARP)+1];
    __shared__ float tempProd[(BLOCKSIZE/HALFWARP)][HALFWARP+PAD];
    __shared__ float ys[INSPECT_BLOCK_c];
    __shared__ int ins_rowInd[2];
    __shared__ int ins_Ind;
        if ((tid <= ((BLOCKSIZE/HALFWARP))) && (myi < numRows))
            rowInd[tid] = rowIndices[myblock+tid];
        __syncthreads();
        if ((ind2Dx < 2) && (ind2Dy == 0) && (myi < numRows)) 
            ins_rowInd[ind2Dx]=ins_rowIndices[bid+ind2Dx];
        __syncthreads();
        t=0;
        unsigned int ind, lb, ub, j;
        if (myi < numRows) {
            lb = rowInd[ind2Dy]+ind2Dx;
            ub = rowInd[ind2Dy+1];
            j = lb;
            if (j<ub) ind = indices[j]; 
            else ind = numCols+INSPECT_BLOCK_c;
        }
        __syncthreads();
        unsigned int ktlb = ins_rowInd[0];
        unsigned int ktub = ins_rowInd[1];
        for (unsigned int kt=ktlb; kt<ktub; kt++) {
            __syncthreads();
            if (tid==0) ins_Ind=ins_indices[kt];
            __syncthreads();
	#if VAR_BLOCK
            unsigned int k = ins_Ind; // In case of var_block, ins_indices 'll have original column index
	#else
            unsigned int k = ins_Ind*INSPECT_BLOCK_c;
	#endif
            if ( (tid < INSPECT_BLOCK_c) && ((k+tid) < numCols) )
                ys[tid] = y[k+tid];
            __syncthreads();
            if (myi < numRows) {
                while ( ind < (k+INSPECT_BLOCK_c) ) {
                    t += val[j] * ys[ind-k];
                    j+=HALFWARP;
                    if (j < ub) ind = indices[j];
                    else { ind = numCols+INSPECT_BLOCK_c; break; }
                }
            }
        }
        tempProd[ind2Dy][ind2Dx] = t;
        __syncthreads();
  #if 0
        if ((ind2Dx == 0) && (myi < numRows)) {
            t=0;
            for (unsigned int k = 0; k<HALFWARP; k++) {
                t += tempProd[ind2Dy][k];
            }
            x[myi] = t;
        }
        __syncthreads();
  #endif
  #if 1
        // Works for HALFWARP=16
        if ((ind2Dx == 0) && (myi < numRows)) {
            t = tempProd[ind2Dy][0] + tempProd[ind2Dy][1] + tempProd[ind2Dy][2] + tempProd[ind2Dy][3] +\
              tempProd[ind2Dy][4] + tempProd[ind2Dy][5] + tempProd[ind2Dy][6] + tempProd[ind2Dy][7] +\
              tempProd[ind2Dy][8] + tempProd[ind2Dy][9] + tempProd[ind2Dy][10]+ tempProd[ind2Dy][11]+\
              tempProd[ind2Dy][12]+ tempProd[ind2Dy][13]+ tempProd[ind2Dy][14]+ tempProd[ind2Dy][15];
            x[myi] = t;
        }
        __syncthreads();
  #endif
  #if 0
        // Works for HALFWARP=16 & 32
        if (!(ind2Dx % 4) && (myi < numRows)) {
            tempProd[ind2Dy][ind2Dx] += tempProd[ind2Dy][ind2Dx+1] + tempProd[ind2Dy][ind2Dx+2] + tempProd[ind2Dy][ind2Dx+3];
        }
        __syncthreads();
        if ((ind2Dx == 0) && (myi < numRows)) {
        #if 1 // for halfwarp 16
            t = tempProd[ind2Dy][0] + tempProd[ind2Dy][4] + tempProd[ind2Dy][8] + tempProd[ind2Dy][12];
        #else // for halfwarp 32
            t = tempProd[ind2Dy][0] + tempProd[ind2Dy][4] + tempProd[ind2Dy][8] + tempProd[ind2Dy][12]+\
               +tempProd[ind2Dy][16] + tempProd[ind2Dy][20] + tempProd[ind2Dy][24] + tempProd[ind2Dy][28];
        #endif
            x[myi] = t;
        }
        __syncthreads();
  #endif
#endif
}

#endif
