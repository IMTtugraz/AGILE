LSQR : LEAST SQUARES SOLVER

This application takes three arguments:
1.) The first argument is the crs matrix filename
    (CRS stands for: compressed row storage)
2.) The second argument is the filename for the right hand side
    vector: y.
3.) The thrid argument is the filename for the resulting x vector

This program solves the equation: Ax = y for x

Example call:
./lsqr.exe A.bin y.bin x.bin


There are some files for testing in the "data" folder.

You can test the c++ lsqr implementation in comparison with
the matlab implementation with: test_lsqr.m
