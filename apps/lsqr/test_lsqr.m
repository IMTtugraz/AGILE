
clear all;
clc;

addpath('../../include/agile/io');


EXECUTABLE = './lsqr.exe';

% file names
CRS_MATRIX_FILE = './data/A_real_64x64_x_75x64.bin';
Y_VECTOR_FILE = './data/y_complex_75x64_nonuniform.bin';
REF_X_VECTOR = './data/x_complex_64x64_uniform.bin';
RES_X_VECTOR = './x_result.bin';

% execute statement
disp(sprintf('execute the following statement:\n%s %s %s %s', EXECUTABLE, CRS_MATRIX_FILE, Y_VECTOR_FILE, RES_X_VECTOR));

% pause
input('press enter when executed ... ');

% run lsqr c++ test
% system(sprintf('%s %s %s %s', EXECUTABLE, MATRIX_FILE, Y_VECTOR_FILE, RES_X_VECTOR));

xref = readbin_vector(REF_X_VECTOR);

% c++ lsqr result
% --------------------------------------------------
x_cpp_res = readbin_vector(RES_X_VECTOR);
diff_cpp = abs(xref - x_cpp_res);
disp(['c++ lsqr max diff: ', num2str(max(diff_cpp))]);


% matlab lsqr result
% --------------------------------------------------
A = readbin_crsmatrix(CRS_MATRIX_FILE);
y = readbin_vector(Y_VECTOR_FILE);
x_matlab_res = lsqr(A, y, 1e-6, 20);
diff_matlab = abs(xref - x_matlab_res);
disp(['matlab lsqr max diff: ', num2str(max(diff_matlab))]);



