SPMVM : SPARSE MATRIX VECTOR MULTIPLICATION

This application takes three arguments:
1.) The first argument is the crs matrix filename
    (CRS stands for: compressed row storage)
2.) The second argument is the vector filename
3.) The thrid argument is the multiplication result vector file name

After loading the matrix and the vector the multiplication
y = Ax 
will be computed.
The result, here described as "y" will be saved to the file given by 
name in the 3rd argument.

Example call:

./spmvm.exe A.bin x.bin y.bin
