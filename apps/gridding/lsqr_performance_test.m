
addpath('src/matlab');
addpath('~/lsqr/');

clear all; clc;

%------------------------------------------------------------------------
% TEST SETTINGS
%------------------------------------------------------------------------
LSQR_TOL = 1e-3;
LSQR_MAXIT = 100;
SHOW_RESULTS_IN_FIGURE = 0;

%------------------------------------------------------------------------
% PERFORMANCE TEST: 4x4
%------------------------------------------------------------------------
% disp('-----------------------------------------------------------------');
% dimX = 4;
% dimY = 4;
% csrMatrixTestFile = './data/csrmatrix_A_16x16.bin';
% kspaceVectorTestFile = './data/vector_kspace_data_16(4x4).bin';
% 
% A = readbin_crsmatrix(csrMatrixTestFile);
% y = readbin_vector(kspaceVectorTestFile);
% 
% tic;
% x = lsqr(A, y, LSQR_TOL, LSQR_MAXIT);
% t = toc;
% 
% if SHOW_RESULTS_IN_FIGURE
%     figure, imshow(abs(ifft2c(reshape(x, dimX, dimY))),[]);
% end
% disp([num2str(dimX),'x', num2str(dimY), ' lsqr time: ', num2str(t), 's']);

%------------------------------------------------------------------------
% PERFORMANCE TEST: 256x256
%------------------------------------------------------------------------
disp('-----------------------------------------------------------------');
dimX = 256;
dimY = 256;
csrMatrixTestFile = './data/csrmatrix_A_76800x65536.bin';
kspaceVectorTestFile = './data/vector_kspace_data_76800(300x256).bin';

A = readbin_crsmatrix(csrMatrixTestFile);
y = readbin_vector(kspaceVectorTestFile);

tic;
x = lsqr(A, y, LSQR_TOL, LSQR_MAXIT);
t = toc;

if SHOW_RESULTS_IN_FIGURE
    figure, imshow(abs(ifft2c(reshape(x, dimX, dimY))),[]);
end
disp([num2str(dimX),'x', num2str(dimY), ' lsqr time: ', num2str(t), 's']);

%------------------------------------------------------------------------
% PERFORMANCE TEST: 512x512
%------------------------------------------------------------------------
disp('-----------------------------------------------------------------');
dimX = 512;
dimY = 512;
csrMatrixTestFile = './data/csrmatrix_A_262144x262144.bin';
kspaceVectorTestFile = './data/vector_kspace_data_262144(512x512).bin';

A = readbin_crsmatrix(csrMatrixTestFile);
y = readbin_vector(kspaceVectorTestFile);

tic;
x = lsqr(A, y, LSQR_TOL, LSQR_MAXIT);
t = toc;

if SHOW_RESULTS_IN_FIGURE
    figure, imshow(abs(ifft2c(reshape(x, dimX, dimY))),[]);
end
disp([num2str(dimX),'x', num2str(dimY), ' lsqr time: ', num2str(t), 's']);

%------------------------------------------------------------------------
% PERFORMANCE TEST: 1024x1024
%------------------------------------------------------------------------
disp('-----------------------------------------------------------------');
dimX = 1024;
dimY = 1024;
csrMatrixTestFile = './data/csrmatrix_A_1048576x1048576.bin';
kspaceVectorTestFile = './data/vector_kspace_data_1048576(1024x1024).bin';

A = readbin_crsmatrix(csrMatrixTestFile);
y = readbin_vector(kspaceVectorTestFile);

tic;
x = lsqr(A, y, LSQR_TOL, LSQR_MAXIT);
t = toc;

if SHOW_RESULTS_IN_FIGURE
    figure, imshow(abs(ifft2c(reshape(x, dimX, dimY))),[]);
end
disp([num2str(dimX),'x', num2str(dimY), ' lsqr time: ', num2str(t), 's']);

return;

%------------------------------------------------------------------------
% PERFORMANCE TEST: 2048x2048
%------------------------------------------------------------------------
% disp('-----------------------------------------------------------------');
% dimX = 2048;
% dimY = 2048;
% csrMatrixTestFile = './data/csrmatrix_A_4194304x4194304.bin';
% kspaceVectorTestFile = './data/vector_kspace_data_4194304(2048x2048).bin';
% 
% A = readbin_crsmatrix(csrMatrixTestFile);
% y = readbin_vector(kspaceVectorTestFile);
% 
% tic;
% x = lsqr(A, y, LSQR_TOL, LSQR_MAXIT);
% t = toc;
% 
% if SHOW_RESULTS_IN_FIGURE
%     figure, imshow(abs(ifft2c(reshape(x, dimX, dimY))),[]);
% end
% disp([num2str(dimX),'x', num2str(dimY), ' lsqr time: ', num2str(t), 's']);
