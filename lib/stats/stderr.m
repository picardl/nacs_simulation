function y = stderr(x,d)
% STDERR    Standard Error
%   For vectors, STDERR(X) is the STD(X)/SQRT(LENGTH(X)). For
%   matrices, STDERR(X) is a row vector containing the standard
%   error value of each column. For N-D arrays, STDERR(X) is the
%   standard error of the elements along the first non-singleton
%   dimension of X.
%
%   STDERR(X,DIM) takes the standard error along the dimension
%   DIM of X.
    if nargin<2
        d=1;
    end
    y = std(shiftdim(x),1,d)./sqrt(size(shiftdim(x),d));
end