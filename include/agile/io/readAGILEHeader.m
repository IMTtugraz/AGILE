function agile_header = readAGILEHeader(fid)
% readAGILEHeader   reads the AGILE header from a file
%
%   agile_header = readAGILEHeader(fid)
%
% Input:
%   fid - The file identifier from which the header is to be read. The
%         file has to be opened before calling this function!
%
% Output:
%   AGILE header as struct with fields
%     .version            - The version as string.
%     .type               - The file type as string.
%     .parameter          - The parameter information as string.
%     .comment            - The file's comment as string.
%     .num_elements       - The number of elements as string.
%     .split_file_type    - Cell array of strings containing the sub-types.
%     .split_num_elements - The number of elements as vector.
%

% $Id: readAGILEHeader.m 461 2011-05-31 14:51:59Z freiberger $

  % prepare the final struct
  agile_header = struct(...
    'version', [], ...
    'type', [], ...
    'split_type', class(0), ...
    'size', [], ...
    'split_size', class(0), ...
    'comment', []);

  % read the version line
  header_line = lower(fgetl(fid));
  if (strfind(header_line, 'version:') ~= 1)
    % error: no version information found on line 1
    error('No file version found on line 1.');
  end
  % cut out the version
  header_line = header_line(9 : end);
  header_line = regexprep(regexprep(header_line, '^ *', ''), ' *$', '');
  % the version has to be 1.0
  if ~strcmp(agile_header.version, '1.0')
    error(['Cannot handle file version "', agile_header.version, '".']);
  end
  agile_header.version = header_line;

  % read the file type
  header_line = lower(fgetl(fid));
  if (strfind(header_line, 'type:') ~= 1)
    % error: no type information found on line 2
    error('No file type found on line 2.');
  end
  header_line = header_line(6 : end);
  header_line = regexprep(regexprep(header_line, '^ *', ''), ' *$', '');
  agile_header.type = header_line;
  % break the file type line into sub-types
  agile_header.split_type = regexp(header_line, ' ', 'split');

  % read the size
  header_line = lower(fgetl(fid));
  if (strfind(header_line, 'size:') ~= 1)
    % error: no size information found on line 3
    error('No size found on line 3.');
  end
  header_line = header_line(6 : end);
  header_line = regexprep(regexprep(header_line, '^ *', ''), ' *$', '');
  agile_header.size = header_line;
  % break the size into sub-sizes
  header_line = regexp(header_line, ' ', 'split');
  agile_header.split_size = zeros(size(header_line));
  for counter = 1 : length(header_line(:))
    agile_header.split_size(counter) = str2double(header_line{counter});
  end

  % read the comment
  header_line = lower(fgetl(fid));
  if (strfind(header_line, 'comment:') ~= 1)
    % error: no comment found on line 4
    error('No comment found on line 4.');
  end
  header_line = header_line(6 : end);
  header_line = regexprep(regexprep(header_line, '^ *', ''), ' *$', '');
  agile_header.comment = header_line;

% End of $Id: readAGILEHeader.m 461 2011-05-31 14:51:59Z freiberger $.
