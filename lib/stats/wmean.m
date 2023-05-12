function [result,err,rchi2] = wmean(m,s,dim)
% WMEAN     Weighted Mean
%   For vectors, WMEAN(X,ERR) is the weighted mean, 
%   SUM(X.*ERR.^-2)./SUM(ERR.^-2). For matrices, WMEAN(X,ERR) is a row 
%   vector containing the weighted mean of each column. For N-D arrays, 
%   WMEAN(X,ERR) is the weighted mean of the elements along the first 
%   non-singleton dimension of X.
%
%   WMEAN(X,ERR,DIM) takes the weighted mean along the dimension
%   DIM of X.
    if nargin<3
        dim = 1;
    end
    if nargin<2
        result = mean(m);
        err = stderr(m);
    else
        m_orig = m;
        
        m = shiftdim(m,dim-1);
        s = shiftdim(s,dim-1);
        
        w = (s+eps).^-2;
        result = sum(w.*m,1)./sum(w,1);

        % Corrected stdev, 160810. This is the unbiased estimate of the
        % sample standard deviation with reliability weights, from
        % https://en.wikipedia.org/wiki/Weighted_arithmetic_mean. A
        % discussion can also be found here:
        % http://mathoverflow.net/questions/11803/unbiased-estimate-of-the-variance-of-a-weighted-mean.
        % Note that the formula derived in that discussion assumes V1=1
        % below.
%         V1 = sum(w,1);
%         V2 = sum(w.^2,1);
%         stdev = sqrt(sum(w.*(m-result).^2,1)./(V1-V2./V1));
%         err = stdev/sqrt(size(m,1));
%         rchi2 = sum(w.*((m-result)./s).^2,1)./(V1-V2./V1);
        
        % conventional error, March 2017
        err = sum(w).^(-1/2);
        rchi2 = sum(((m-result)./s).^2)./(size(m,1)-1);
        
        if numel(m)==1
            result = m;
            err = s;
        end
        
        result = shiftdim(result,length(size(m_orig))-dim+1);
        err = shiftdim(err,length(size(m_orig))-dim+1);
        rchi2 = shiftdim(rchi2,length(size(m_orig))-dim+1);
    end
end