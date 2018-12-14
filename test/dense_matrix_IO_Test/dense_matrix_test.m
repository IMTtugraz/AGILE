% Skript: Test Matlab/CUDA Binary File I/O
%
% NOTE: Calling fft_test.exe might not work due to problems with libstdc++:
% In this case, you have to run fft_test.exe from the command line:
% ./dense_matrix_IO_test.exe testfile1.bin testfile2.bin
% 
% Last Change: 26.3.2009, 14:30, by Florian

clear all, close all, clc;

addpath('/home/knoll/svn/imtcuda/include/imtcuda/io'); % Change this

%% Control Structures
testComplex = 1;

%% generate test matrix

if testComplex
    AReal = [1:5; 6:10; 11:15; 16:20];
    AImag = [i*(21:25); i*(26:30); i*(31:35); i*(36:40)];
    A = single(AReal+AImag)
else
    A = single([1:5; 6:10; 11:15; 16:20])
end

[numRows,numColums] = size(A);
numBytesPerEntry = 4;

%% Start with File I/O
% write image data to a file
writeMatlab2Bin(A, 'testfile1.bin');

% Just to check: Read in just right now
A2 = readBin2Matlab('testfile1.bin', numRows, numColums, numBytesPerEntry)

% Run file that reads binary file in AGILE and writes result back
% If it does not work, run
% ./dense_matrix_IO_test.exe testfile1.bin testfile2.bin
% on the console
system('./dense_matrix_IO_test.exe testfile1.bin testfile2.bin');

% Read in result that was written by imtcuda
A3 = readBin2Matlab('testfile2.bin', numRows, numColums, numBytesPerEntry)
