function kshow(kspace_data, gamma, stretch)

%  Plot scaled absolute value of kspace, logarithmic scale
%========================================================================
%  KSHOW(KSPACE_DATA)
%========================================================================
%
% This syntax is used a lot, so it is saved as a function to avoid typing
% overhead
if nargin<3
    stretch=[];
end

if nargin<2
    gamma=0;
end

%imshow(log(1+abs(kspace_data)),[]) % Default

imshow(brighten(log(1+abs(kspace_data)),gamma),stretch)