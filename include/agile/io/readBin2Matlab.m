function matrix = readBin2Matlab(filename, numRows, numColumns, numBytesPerEntry)
% MATRIX = READBIN2MATLAB(FILENAME, NUMROWS, NUMCOLUMS, NUMBYTESPERENTRY)
% Reads a bin file from MTCUDA in matlab
%
% INPUT:
% filename: string, name of the bin file
% numRows 
% numColums
% numBytesPerEntry
%
% OUTPUT:
% matrix: any 2D matrix, Single precision, real or complex
%
% Example: testMatrix = readBin2Matlab('testfile1.bin', numRows,
%                       numColums, numBytesPerEntry);
%
% See also: writeMatlab2bin.m
%
% Last Change: 26.3.2009, 14:00
% By: Florian

%% open file
fid = fopen(filename, 'rb');
if fid == -1
  error(['Could not open file ', filename '']);
end

%% Check file information
matrix_info.rows = fread(fid, 1, 'uint32');
if matrix_info.rows ~= numRows
  error('Dimension mismatch');
end
matrix_info.columns = fread(fid, 1, 'uint32');
if matrix_info.columns ~= numColumns
  error('Dimension mismatch');
end
matrix_info.num_byte_per_entry = fread(fid, 1, 'uint32');
if matrix_info.num_byte_per_entry ~= numBytesPerEntry
  error('Type mismatch');
end

%% Determine if file contains complex or real content
matrix_info.is_complex = fread(fid, 1, 'uint32');
if matrix_info.is_complex
  temp = single(fread(fid, matrix_info.rows * matrix_info.columns * 2, 'float32'));
  matrix = temp(1 : 2 : end) + i * temp(2 : 2 : end);
else
  matrix = single(fread(fid, matrix_info.rows * matrix_info.columns, 'float32'));
end

%% Arrange in the right shape: Note that matlab permutes rows and columns
%% in comparison to C. This is why the reshape command has to be called that way 
fclose(fid);
matrix = reshape(matrix,numColumns,numRows).';
