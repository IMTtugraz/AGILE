#!/bin/bash

make

# clear old log files
#------------------------------------------------------------------------------
rm -f test_results_*.log

# test with different source images and increasing number of interp. points
#------------------------------------------------------------------------------
for (( i = 512;  i <= 2048;  i=i*2 ))
do
  for (( j=10000; j<=1000000; j=j*10 ))
  do
    for (( x=0, v=j/2; x<2; x++, v=j ))
    do
      ./interp2d.exe $i $i $v | tee -a test_results_${i}x${i}.log
    done
  done
done


# tests with non-square source images
#------------------------------------------------------------------------------

# wide images (width > texture cache size)
#------------------------------------------------------------------------------
NUM_COLS=16384
NUM_ROWS=512
for (( j=10000; j<=1000000; j=j*10 ))
do
  for (( x=0, v=j/2; x<2; x++, v=j ))
  do
    ./interp2d.exe $NUM_COLS $NUM_ROWS $v | tee -a test_results_${NUM_COLS}x${NUM_ROWS}.log
  done
done

NUM_COLS=8192
NUM_ROWS=1024
for (( j=10000; j<=1000000; j=j*10 ))
do
  for (( x=0, v=j/2; x<2; x++, v=j ))
  do
    ./interp2d.exe $NUM_COLS $NUM_ROWS $v | tee -a test_results_${NUM_COLS}x${NUM_ROWS}.log
  done
done


# tall images (width > texture cache size)
#------------------------------------------------------------------------------
NUM_COLS=512
NUM_ROWS=16384
for (( j=10000; j<=1000000; j=j*10 ))
do
  for (( x=0, v=j/2; x<2; x++, v=j ))
  do
    ./interp2d.exe $NUM_COLS $NUM_ROWS $v | tee -a test_results_${NUM_COLS}x${NUM_ROWS}.log
  done
done

NUM_COLS=1024
NUM_ROWS=8192
for (( j=10000; j<=1000000; j=j*10 ))
do
  for (( x=0, v=j/2; x<2; x++, v=j ))
  do
    ./interp2d.exe $NUM_COLS $NUM_ROWS $v | tee -a test_results_${NUM_COLS}x${NUM_ROWS}.log
  done
done


