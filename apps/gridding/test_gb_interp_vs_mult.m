%--------------------------------------------------------------------------
% BACKWARD GRIDDING TEST - Interpolation vs. Multiplication
%--------------------------------------------------------------------------

% Config
%--------------------------------------------------------------------------
TRANSFOMATION_MATRIX_FILE = 'data/csrmatrix_AT_16x16.bin';
KSPACE_POS_VECTOR_FILE = 'data/vector_kspace_positions_16(4x4).bin';
GSPACE_SOURCE_VECTOR_FILE = 'data/vector_gspace_data_16(4x4).bin';

INTERPOLATION_RESULT_FILE = 'res_interp.bin';
MULTIPLICATION_RESULT_FILE = 'res_mult.bin';

% Tests
%--------------------------------------------------------------------------
clc;

% read files
A = readbin_crsmatrix(TRANSFOMATION_MATRIX_FILE);
pos = reshape(readbin_vector(KSPACE_POS_VECTOR_FILE), 4, 4) + 1;
ri = reshape(readbin_vector(INTERPOLATION_RESULT_FILE), 4,4);
rm = reshape(readbin_vector(MULTIPLICATION_RESULT_FILE), 4,4);
src = reshape(readbin_vector(GSPACE_SOURCE_VECTOR_FILE), 4, 4)';

% src = complex(1:16,0);
% src = reshape(src, 4, 4)'
% interp2(src, 2, 1)
% return;

clear tmp;
tmp(:,1) = pos(:);
tmp(:,2) = rm(:);
tmp(:,3) = ri(:);
tmp(:,4) = abs(rm(:)-ri(:)) / max(abs(ri(:)));
tmp(:,5) = A*src(:);
tmp(:,6) = interp2(src, real(pos(:)), imag(pos(:)));

tmp
