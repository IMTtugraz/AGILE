#!/bin/bash

#########################################
# ATTENTION                             #
# ibm spmv does not work with cuda2.3   #
# (no compile errors but segfault)      #
#########################################

make
make -C spmv/

IMT_EXECUTABLE="./comp_spmv.exe"
IBM_EXECUTABLE="./spmv/SpMV"

LOG_FILE="./spmv.log"

MATRICES=(
  "rail4284.mtx" # too big for 9600M GT
  "rdist1.mtx"
  "orani678.mtx"
  "rim.mtx"
  "rma10.mtx"
  "para-7.mtx"
  "e40r0100.mtx"
  "olafu.mtx"
  "bcsstk35.mtx"
  "venkat01.mtx"
  "nasasrb.mtx"
  "ex11.mtx"
  "raefsky3.mtx"
  "interp512.mtx"
  "interp1024.mtx"
)



#-----------------------------------------------------------------------------

rm -f $LOG_FILE;
echo -ne "# Matrixname\tMatrix rows\tMatrix columns\tIMT GPU (ms)\tIMT CPU (ms)\tIMT GFLOPS\tIMT Success\tIBM GPU (ms)\tIBM CPU (ms)\tIBM GFLOPS\tIBM Success" | tee -a $LOG_FILE
for mat in ${MATRICES[@]}
do
  echo -ne "\n#----------------------------------------> $mat\n" | tee -a $LOG_FILE
  echo -ne "$mat\t" | tee -a $LOG_FILE
  $IMT_EXECUTABLE "matrices/$mat" 2>&1 | tee -a $LOG_FILE
  echo -ne "\t" | tee -a $LOG_FILE
  $IBM_EXECUTABLE "matrices/$mat" 2>&1 | tee -a $LOG_FILE
done
echo -ne "\n" | tee -a $LOG_FILE
