% Clean all test data produced by "build_test_data.m"
% Gerald Buchgraber <gerald.buchgraber@student.tugraz.at>
%==========================================================================

% gridding data target folder
targetFolder = './data';

delete([targetFolder '/vector_kspace_positions_*.bin']);
delete([targetFolder '/csrmatrix_A_*.bin']);
delete([targetFolder '/csrmatrix_AT_*.bin']);