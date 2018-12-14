function result = readAGILEVariable(file_name)
% readAGILEVariable   read an AGILE variable from a file
%
%   result = readAGILEVariable(file_name)
%
% Input:
%   file_name - the name of the file to read
%
% Output:
%   the variable read from the file
%

% $Id: readAGILEVariable.m 460 2011-05-31 14:51:17Z freiberger $

  fid = fopen(file_name);
  if (fid == -1)
    error(['Could not open file "', file_name, '".']);
  end
  header = readAGILEHeader(fid);

  % handle different file types
  switch(header.split_type{1})
    case 'matrix'
      % branch depending on the matrix type
      switch(header.split_type{2})
        case 'dense'
          result = readDenseMatrix(fid, header);
        otherwise
          error(['Cannot read matrix type "', header.split_type{2}, '".']);
      end

    case 'vector'
      % branch depending on the vector type
      switch(header.split_type{2})
        case 'dense'
          result = readDenseVector(fid, header);
        otherwise
          error(['Cannot read vector type "', header.split_type{2}, '".']);
      end

    otherwise
      error(['Cannot handle file type "', header.type, '".']);
  end
  fclose(fid);
end


% read a dense matrix
function result = readDenseMatrix(fid, header)
  is_complex = strcmp(header.split_type{3}, 'complex');
  num_rows = header.split_size{1};
  num_columns = header.split_size{2};
  if ~is_complex
    temp = fread(fid, [num_columns, num_rows], ['*',header.split_type{4}]);
    result = temp.';
  else
    temp = fread(fid, [num_columns * 2, num_rows], ...
      ['*', header.split_type{4}]);
    result = (temp(1:2:end, :) + i*temp(2:2:end)).';
  end
end


% read a dense vector
function result = readDenseVector(fid, header)
  is_complex = strcmp(header.split_type{3}, 'complex');
  vector_size = header.split_size{1};
  if ~is_complex
    result = fread(fid, vector_size, ['*', header.split_type{4}]);
  else
    temp = fread(fid, 2 * vector_size, ['*', header.split_type{4}]);
    result = temp(1:2:end) + i*temp(2:2:end);
  end
end

% End of $Id: readAGILEVariable.m 460 2011-05-31 14:51:17Z freiberger $.
