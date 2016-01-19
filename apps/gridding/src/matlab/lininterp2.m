function A = lininterp2(arg1,arg2,arg3,arg4)
%LININTERP2 2-D bilinear data interpolation.
%   A = LININTERP2(X,Y,XI,YI) uses bilinear interpolation to
%   build the transformation matrix A used to transform data corresponding
%   to the points given in X and Y to interpolated data corresponding to
%   the points given in XI and YI.
%
%   Matrices X and Y specify the points at which the transformation input
%   data is located.
%   X and Y can also be vectors specifying the abscissae as for MESHGRID.
%   In both cases, X and Y must be equally spaced and monotonic.
%
%   If XI and YI are vectors, LINEAR returns vector ZI containing
%   the interpolated values at the corresponding points (XI,YI).
%


if nargin==2 % lininterp2(x,y), expand data
    ncols = length(arg1);
    nrows = length(arg2);
    s = 1:.5:ncols; lengths = length(s);
    t = (1:.5:nrows)'; lengtht = length(t);
    s = repmat(s,lengtht,1);
    t = repmat(t,1,lengths);

elseif nargin==3 % lininterp2(x,y,n), expand data n times
    [nrows,ncols] = size(arg1);
    ntimes = floor(arg2);
    s = 1:1/(2^ntimes):ncols; lengths = length(s);
    t = (1:1/(2^ntimes):nrows)'; lengtht = length(t);
    s = repmat(s,lengtht,1);
    t = repmat(t,1,lengths);
    
elseif nargin==4 % lininterp2(x,y,xi,yi), X and Y specified.
    if (~isequal(size(arg1),size(arg2)))
        error('MATLAB:interp2:linear:XYZLengthMismatch',...
            'The lengths of the X and Y vectors must match Z.');
    end
    
    [nrows, ncols] = size(arg1);
    if nrows < 2 || ncols < 2
        error('lininterp2:sizeXY','X and Y must be at least 2-by-2.');
    end

    s = 1 + (arg3-arg1(1))/(arg1(end)-arg1(1))*(ncols-1);
    t = 1 + (arg4-arg2(1))/(arg2(end)-arg2(1))*(nrows-1);

elseif nargin < 2 || nargin > 4
    error('lininterp2:nargin','Wrong number of input arguments.');
end


if nrows < 2 || ncols < 2
    error('lininterp2:sizeXY','X and Y must be at least 2-by-2.');
end

if ~isequal(size(s),size(t))
    error('lininterp2:XIandYISizeMismatch',...
        'XI and YI must be the same size.');
end

% Check for out of range values of s and set to 1
sout = find((s<1)|(s>ncols));
if ~isempty(sout), s(sout) = 1; end

% Check for out of range values of t and set to 1
tout = find((t<1)|(t>nrows));
if ~isempty(tout), t(tout) = 1; end

% Matrix element indexing
ndx = floor(t)+floor(s-1)*nrows;

% Compute intepolation parameters, check for boundary value.
if isempty(s), d = s; else d = find(s==ncols); end
s(:) = (s - floor(s));
if ~isempty(d), s(d) = s(d)+1; ndx(d) = ndx(d)-nrows; end

% Compute intepolation parameters, check for boundary value.
if isempty(t), d = t; else d = find(t==nrows); end
t(:) = (t - floor(t));
if ~isempty(d), t(d) = t(d)+1; ndx(d) = ndx(d)-1; end

% Now interpolate.
ndx = ndx(:);
sizeS = size(s);
sizeX = [ncols nrows];

s = s(:);
t = t(:);
onemt = 1-t;

idx = 1:length(ndx);
rows = repmat(idx,1,4)';
cols = [ndx(idx); ndx(idx)+1; ndx(idx)+nrows; ndx(idx)+(nrows+1)];

vals = [onemt(idx).*(1-s(idx)); t(idx).*(1-s(idx)); onemt(idx).*s(idx); t(idx).*s(idx)];

A = sparse(rows, cols, vals, sizeS(1)*sizeS(2), sizeX(1)*sizeX(2));
