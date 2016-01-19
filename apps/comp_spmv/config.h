/*
 * IBM Sparse Matrix-Vector Multiplication Toolkit for Graphics Processing Units
 *  (c) Copyright IBM Corp. 2008, 2009.  All Rights Reserved.
 */ 

#ifndef __CONFIG_H__
#define __CONFIG_H__

#define INSPECT 0
#define INSPECT_BLOCK_c (BLOCKSIZE)  
#define VAR_BLOCK 1
#define VAR_COLUMN 32 

#define INSPECT_INPUT 0
#define INSPECT_INPUT_MAX 512

#define NUMTHREADS 512
#define BLOCKSIZE 512
#define HALFWARP 16

#define VERIFY 1
#define DEBUG_R 0
#define EXEC_CPU 1
#define NUM_ITER 10

#define BCSR 0
#define BCSR_r 8
#define BCSR_c 8
#define BCSR_c8 1
#define BCSR_c4 0
#define BCSR_c2 0
#define BCSR_r2 0 

#define PADDED_CSR 1

#define C_GLOBAL_OPT 0
#define NEW_GLOBAL_OPT 0 
#define GLOBAL_OPT 1
#define GLOBAL_SHARED_OPT 0
#define GLOBAL 0
#define SHARED_RI 0
#define PAD 1

#define CACHE 1

#define TIMER 1

#endif /* __CONFIG_H__ */
