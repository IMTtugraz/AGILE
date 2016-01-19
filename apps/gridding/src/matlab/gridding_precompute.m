function [ ] = gridding_precompute( folder, gResX, gResY, kNumTraj, kNumSamplesPerTraj )
%GRIDDING_PRECOMPUTE - Precompute necessary gridding data files
%   The files that will be precomputed are:
%       - the transformation matrix A
%       - the adjoint transformation matrix A'
%       - and the position vector of the kspace sample points
%
%   The transformation matrix A describes the bininear interpolation.

if ~ischar(folder)
    error('The first parameter (target folder) must be a string.');
elseif ~exist(folder, 'dir')
    error(['Folder: "' folder '" does not exist.']);
end

if nargin < 2 || nargin == 4
    error('Wrong number of input arguments.');
elseif nargin == 2
    gResY = gResX;
    kNumTraj = gResX;
    kNumSamplesPerTraj = gResX;
elseif nargin == 3
    kNumTraj = gResY;
    kNumSamplesPerTraj = gResY;
    gResY = gResX;
end

% gridded space (positions)
%========================================================================
% this was the old algorithm to form grid coords
% rhoX = linspace(-0.5, 0.5, gResX);
% rhoY = linspace(0.5, -0.5, gResY)';
% 
% x = zeros(gResY, gResX);
% for ii = 1:gResY
%     x(ii,:) = rhoX + j*rhoY(ii);
% end
% 
% scale and move x to be between (0,0) and (gResX, gResY)
% x = complex((real(x)+0.5) * (gResX-1), abs(imag(x)-0.5) * (gResY-1));

[xt yt] = meshgrid(0:1:(gResX-1), 0:1:(gResY-1));
x = complex(xt, yt);

% radial sampling (positions)
%========================================================================
theta = 0:pi/kNumTraj:pi-pi/kNumTraj;
rho = linspace(-0.5, 0.5, kNumSamplesPerTraj)';

b = zeros(kNumSamplesPerTraj, kNumTraj);
for ii = 1:length(theta)
    b(:,ii) = rho*exp(-j*theta(ii));
end
% scale and move b to be between (0,0) and (gResX, gResY)
b = complex((real(b)+0.5) * (gResX-1), abs(imag(b)-0.5) * (gResY-1));


% quantization of texture coords as gpu does
b = round(b * 256) / 256;

% transformation matrix (describing the bilinear interpolation)
A = lininterp2(real(x), imag(x), real(b), imag(b));

% imaginary part of x and b is multiplicated by -1 for better visualization
%   => GPU addresses pixel with 0,0 in upper left corner
% x = complex(real(x), imag(x) * -1);
% b = complex(real(b), imag(b) * -1);
% bMat = reshape(b, kNumSamplesPerTraj, kNumTraj);
% xMat = reshape(x, gResY, gResX);
% figure
% subplot(2,2,1), plot(bMat, 'r-'), hold, plot(real(bMat), imag(bMat), 'rx', real(b(1)), imag(b(1)), 'ko', real(b(2)), imag(b(2)), 'ks', real(bMat(1,2)), imag(bMat(1,2)), 'k^'), title('Non Uniform Grid');
% subplot(2,2,2), plot(real(xMat), imag(xMat), 'bx', real(x(1)), imag(x(1)), 'ko', real(x(2)), imag(x(2)), 'ks', real(bMat), imag(bMat), 'rx'), title('Uniform Grid (blue)');
% subplot(2,2,4), spy(A), title('Transformation Matrix (spy)');


% save kspace position vector
writebin_vector(b(:), [folder sprintf('/vector_kspace_positions_%d(%dx%d).bin', kNumTraj*kNumSamplesPerTraj, kNumTraj, kNumSamplesPerTraj)]);



% save transformation matrix
writebin_crsmatrix(A, [folder sprintf('/csrmatrix_A_%dx%d.bin', kNumTraj*kNumSamplesPerTraj, gResX*gResY)]);

% save adjoiont transformation matrix
writebin_crsmatrix(A', [folder sprintf('/csrmatrix_AT_%dx%d.bin', gResX*gResY, kNumTraj*kNumSamplesPerTraj)]);

end
