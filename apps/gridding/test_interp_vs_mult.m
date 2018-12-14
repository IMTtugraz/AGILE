clear all; clc;

dimX = 4;
dimY = 4;

A = readbin_crsmatrix('data/csrmatrix_A_16x16.bin');
x = readbin_vector('data/vector_gspace_data_16(4x4).bin');
p = readbin_vector('data/vector_kspace_positions_16(4x4).bin');
% fix positions for matlab (starting from index one)
p = p + 1 + i;
xd = x;
yi = readbin_vector('res_interp.bin');
ym = readbin_vector('res_mult.bin');

% convert to single precision
yi = single(yi);
ym = single(ym);
A = single(full(A));
x = single(x);
p = single(p);

% compute gold solution
ym_gold = A*x;
yi_gold = test_interp2(reshape(x, dimY, dimX), real(p), imag(p));

% cuda texture fetch (2D)
% tex( x, y ) = (1 − α )(1 − β )T [i, j ] + α (1 − β )T [i + 1, j ] + (1 − α ) β T [i, j + 1] + αβ T [i + 1, j + 1]
%
% α, β  ... 9-bit fixed point format with 8 bits of fractional value
% 


% tpIdx = 6;
% tri = test_interp2(reshape(x, dimY, dimX), real(p(tpIdx)), imag(p(tpIdx)))
% trm = A(tpIdx,:)*x
% tri-trm


% yi_gold - ym_gold
% disp([ '' num2str(norm(yi_gold - ym_gold)) ]);
