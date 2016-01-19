% Precompute:
%   - bilinear interpolation transformation matrix (and adjoint)
%   - position vector for radial samples (required for gpu interpolation)
%
% Gerald Buchgraber <gerald.buchgraber@student.tugraz.at>
%==========================================================================

addpath('src/matlab');

% gridding data target folder
targetFolder = './data';

% gridding_precompute(targetFolder, 4);

gridding_precompute(targetFolder, 256, 256, 300, 256);

gridding_precompute(targetFolder, 512);

gridding_precompute(targetFolder, 1024);

%gridding_precompute(targetFolder, 2048);


