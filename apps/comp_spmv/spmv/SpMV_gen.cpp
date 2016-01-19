/*
 * IBM Sparse Matrix-Vector Multiplication Toolkit for Graphics Processing Units
 * (c) Copyright IBM Corp. 2008, 2009.  All Rights Reserved.
 */ 

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "SpMV.h"
#include "../config.h"

int cmpRow(const void *e1, const void *e2) {
    if (((NZEntry *)e1)->rowNum < ((NZEntry *)e2)->rowNum) return -1;
    if (((NZEntry *)e1)->rowNum > ((NZEntry *)e2)->rowNum) return 1;
    if (((NZEntry *)e1)->colNum < ((NZEntry *)e2)->colNum) return -1;
    if (((NZEntry *)e1)->colNum > ((NZEntry *)e2)->colNum) return -1;
    return 0;
}

int cmpCol(const void *e1, const void *e2) {
    if (((NZEntry *)e1)->colNum < ((NZEntry *)e1)->colNum) return -1;
    if (((NZEntry *)e1)->colNum > ((NZEntry *)e2)->colNum) return -1;
    if (((NZEntry *)e1)->rowNum < ((NZEntry *)e2)->rowNum) return -1;
    if (((NZEntry *)e1)->rowNum > ((NZEntry *)e2)->rowNum) return 1;
    return 0;
}

void readInputVector(float *y, const char *filename, int numCols)
{
    FILE *f;
    f = fopen(filename,"r");
    if (!f) {
        fprintf(stderr,"Cannot open file: %s\n",filename);
        exit(-1);
    }

    int i=0;
    char line[256];
    while ( ((fgets(line, 256, f)) != NULL) && (i<numCols) ) {
	sscanf(line,"%f", &(y[i]));
	i++;
    }
    if (i<numCols) exit(-1);
    fclose(f);
}

void writeOutputVector(float *x, const char *filename, int numRows)
{
    FILE *f;
    f = fopen(filename,"w");
    if (!f) {
        fprintf(stderr,"Cannot open file: %s\n",filename);
        exit(-1);
    }

    int i;
    for (i = 0; i < numRows; i++)
        fprintf(f,"%f\n", x[i]);

    fclose(f);
}

void readSparseMatrix(SpMatrix *m, const char *filename, int format)
{
    FILE *f;
    f = fopen(filename,"r");
    if (!f) {
        fprintf(stderr,"Cannot open file: %s\n",filename);
        exit(-1);
    }

    char line[256];
    while ( (fgets(line, 256, f)) != NULL) {
	if (line[0] != '%') break;  
    }

    if ( (sscanf(line,"%d %d %d", &(m->numRows), &(m->numCols), &(m->numNZEntries))) != 3) exit(-1);
    m->nzentries = (NZEntry *) malloc(sizeof(NZEntry) * (m->numNZEntries));
    m->rowPtrs = (unsigned int *) malloc(sizeof(unsigned int) * (m->numRows));
    m->colPtrs = (unsigned int *) malloc(sizeof(unsigned int) * (m->numCols));

    NZEntry e;
    int i;
    for (i = 0; (unsigned)i < m->numNZEntries; i++)
    {
        fscanf(f,"%d %d %f\n", &(e.rowNum), &(e.colNum), &(e.val));
	e.rowNum--; e.colNum--;
        (m->nzentries)[i]= e;
    }

    // sort into row-major order or column major order based on the format
    if (format == 0) // ROW-MAJOR
        qsort(m->nzentries, m->numNZEntries, sizeof(NZEntry), cmpRow);
    else if (format == 1) // COLUMN-MAJOR
        qsort(m->nzentries, m->numNZEntries, sizeof(NZEntry), cmpCol);
    else { }

    if (format == 0) {
    	// set index of first elt in each row
    	// relies on at least one item in each row
     	m->rowPtrs[0]=0;
    	unsigned int row, prevrow=0;
        for (i = 1; i < m->numNZEntries; i++) {
            row = (m->nzentries)[i].rowNum;
            if (row != prevrow) {
		prevrow = row;
		m->rowPtrs[prevrow]=i;
	    }
        }
    }

    if (format == 1) {
    	// set index of first elt in each col
    	// relies on at least one item in each col
     	m->colPtrs[0]=0;
    	unsigned int col, prevcol=0;
        for (i = 1; i < m->numNZEntries; i++) {
            col = (m->nzentries)[i].colNum;
            if (col != prevcol) {
		prevcol = col;
		m->colPtrs[prevcol]=i;
	    }
        }
    }
    fclose(f);
}

void genCSRFormat(SpMatrix * m, float *val, unsigned int *rowIndices, unsigned int *indices)
{
    unsigned int numRows = m->numRows;
    for (unsigned int i = 0; i < m->numNZEntries; i++) {
        val[i] = (m->nzentries)[i].val;
        indices[i] = (m->nzentries)[i].colNum;
    }

    for (unsigned int i = 0; i < numRows; i++) {
	rowIndices[i] = m->rowPtrs[i];
    }
    rowIndices[numRows] = m->numNZEntries;
}

void genCSCFormat(SpMatrix * m, float *val, unsigned int *colIndices, unsigned int *indices)
{
    unsigned int numCols = m->numCols;
    for (unsigned int i = 0; i < m->numNZEntries; i++) {
        val[i] = (m->nzentries)[i].val;
        indices[i] = (m->nzentries)[i].rowNum;
    }

    for (unsigned int i = 0; i < numCols; i++) {
        colIndices[i] = m->colPtrs[i];
    }
    colIndices[numCols] = m->numNZEntries;
}

void genBCSRFormat(SpMatrix *m, float **val, unsigned int **rowIndices, unsigned int **indices,
                   unsigned int *numblocks, unsigned int bsx, unsigned int bsy)
{
    unsigned int nblocks=0;
    unsigned int nrows = m->numRows;
    unsigned int ncols = m->numCols;
    unsigned int *lb = (unsigned int *)malloc(sizeof(int)*bsx);
    unsigned int *ub = (unsigned int *)malloc(sizeof(int)*bsx);
    unsigned int *bptr = (unsigned int *)malloc(sizeof(int)*bsx);
    unsigned int *colFlag = (unsigned int *)malloc(sizeof(int)*(int)ceild(ncols,bsy));
    unsigned int *nblocksRow = (unsigned int *)malloc(sizeof(int)*(int)ceild(nrows,bsx));
    unsigned int **indRows = (unsigned int **)malloc(sizeof(int*)*(int)ceild(nrows,bsx));

    *rowIndices = (unsigned int *)malloc(sizeof(int)*((int)ceild(nrows,bsx)+1));
   

    // for each block of row 
    unsigned int it, iti;
    for (it = 0, iti = 0; it < nrows; it += bsx, iti++) {
	// start of a row block
	(*rowIndices)[iti]=nblocks;
        nblocksRow[iti]=0;
        for (unsigned int i = it; i < min(it+bsx,nrows); i++) {
            lb[i-it] = (m->rowPtrs)[i];
            if (i==(nrows-1)) ub[i-it] = m->numNZEntries-1;
            else ub[i-it] = (m->rowPtrs)[i+1]-1;
	    bptr[i-it] = lb[i-it];
	}
        // for each block of column within a row block
	for (unsigned int jt = 0, jti = 0; jt < ncols; jt += bsy, jti++) {
	    colFlag[jti]=0;
  	    unsigned int blockStart = nblocks;
	    for (unsigned int i = it; i < min(it+bsx,nrows); i++) {
		unsigned int j = bptr[i-it];
		for (; j <= ub[i-it]; j++) {
		    unsigned cInd = (m->nzentries)[j].colNum;
		    if (cInd >= jt+bsy)  
			break;
		    if (blockStart == nblocks) {
			nblocks++;
			nblocksRow[iti]++;
			colFlag[jti]=1;
		    }
		}
		bptr[i-it] = j; 
	    }
	}
        indRows[iti] = (unsigned int *)malloc(sizeof(int)*nblocksRow[iti]);
        for (unsigned int k=0, indRowk=0; k < ceild(ncols,bsy); k++) {
            if (colFlag[k]) {
                indRows[iti][indRowk]=k;
                indRowk++;
            }
        }
	
    }
    (*rowIndices)[iti]=nblocks;

    *numblocks = nblocks;
    *indices = (unsigned int *)malloc(sizeof(int)*nblocks);
    *val = (float *)malloc(sizeof(float)*(nblocks*bsx*bsy));

    // Merge all indRows to generate indices
    nblocks=0;
    for (unsigned int k=0; k < ceild(nrows,bsx); k++) {
	for (unsigned int l=0; l<nblocksRow[k]; l++) {
	    (*indices)[nblocks]=indRows[k][l];
	    nblocks++;
	}
    }
    for (unsigned int k=0; k < ceild(nrows,bsx); k++)
	free(indRows[k]);
    free(nblocksRow);
    free(indRows);

    // One more loop to fill in val 
    nblocks=0;
    for (it = 0, iti = 0; it < nrows; it += bsx, iti++) {
	unsigned int lbb = (*rowIndices)[iti];
	unsigned int ubb = (*rowIndices)[iti+1]-1;
        for (unsigned int i = it; i < min(it+bsx,nrows); i++) {
            lb[i-it] = (m->rowPtrs)[i];
            if (i==(nrows-1)) ub[i-it] = m->numNZEntries-1;
            else ub[i-it] = (m->rowPtrs)[i+1]-1;
            bptr[i-it] = lb[i-it];
	}
        for (unsigned int jb = lbb; jb <= ubb; jb++) {
	    unsigned int jti = (*indices)[jb];
            for (unsigned int i = it; i < min(it+bsx,nrows); i++) {
		for (unsigned int k = 0; k < bsy; k++)
		    (*val)[nblocks*bsx*bsy+(i-it)*bsy+k]=0;
                unsigned int j = bptr[i-it];
                for (; j <= ub[i-it]; j++) {
                    unsigned cInd = (m->nzentries)[j].colNum;
                    if (cInd >= ((jti*bsy)+bsy)) break;
		    else (*val)[nblocks*bsx*bsy+(i-it)*bsy+(cInd-(jti*bsy))]=(m->nzentries)[j].val;
		}
		bptr[i-it] = j;
            }
	    nblocks++;
        }
    }

    free(lb);
    free(ub);
    free(bptr);
    free(colFlag);
}

void genPaddedCSRFormat(SpMatrix * m, float **val, unsigned int **rowIndices, unsigned int **indices)
{
    unsigned int numRows = m->numRows;
    (*rowIndices) = (unsigned int *)malloc(sizeof(int)*(numRows+1));
    unsigned int prevRowNum = -1;
    unsigned int padNZ=0, nnzRow, padNZCnt=0;

    for (unsigned int i = 0; i < numRows-1; i++) {
        nnzRow = m->rowPtrs[i+1]-m->rowPtrs[i];
	if (nnzRow%HALFWARP) {
	    padNZ += nnzRow + HALFWARP-(nnzRow%HALFWARP); 
	    //printf("Row:%d --- OrigNZ: %d -- PadNZ: %d\n",i,nnzRow,(nnzRow + HALFWARP-(nnzRow%HALFWARP)));
	}
	else {
	    padNZ += nnzRow;
	    //printf("Row:%d --- OrigNZ: %d -- PadNZ: %d\n",i,nnzRow,nnzRow);
	}
    }
    nnzRow = m->numNZEntries - m->rowPtrs[numRows-1];
    if (nnzRow%HALFWARP) {
        padNZ += nnzRow + HALFWARP-(nnzRow%HALFWARP);
        //printf("Row:%d --- OrigNZ: %d -- PadNZ: %d\n",numRows-1,nnzRow,(nnzRow + HALFWARP-(nnzRow%HALFWARP)));
    }
    else {
        padNZ += nnzRow;
        //printf("Row:%d --- OrigNZ: %d -- PadNZ: %d\n",numRows-1,nnzRow,nnzRow);
    }

    (*val) = (float *)malloc(sizeof(float)*padNZ);
    (*indices) = (unsigned int *)malloc(sizeof(int)*padNZ);

    for (unsigned int i = 0; i < m->numNZEntries; i++) {
	unsigned int currRowNum = (m->nzentries)[i].rowNum;
	if (currRowNum != prevRowNum) { //start of a row
	    if (currRowNum && (padNZCnt%HALFWARP)) { //Not first row and padNZCnt not a multiple of HALFWARP
 		unsigned int fillCount = HALFWARP-(padNZCnt%HALFWARP);
		for (unsigned int j=0;j<fillCount;j++) {
		    (*val)[padNZCnt]=0;
		    (*indices)[padNZCnt]=0;
		    padNZCnt++;
		}
	    }
	    (*rowIndices)[currRowNum]=padNZCnt;
	    prevRowNum = currRowNum;
	}
        (*val)[padNZCnt] = (m->nzentries)[i].val;
        (*indices)[padNZCnt] = (m->nzentries)[i].colNum;
	padNZCnt++;
    }

    (*rowIndices)[numRows] = padNZCnt;
}


