function [y_out,err_out] = prop_errors(f,y,err)
% PROP_ERRORS     Propagates errors through a function
    y_out = mat2cell(y,size(y,1),ones(1,size(y,2)));
    y_out = f(y_out{:});
    
    err_out = zeros(size(y_out));
    for i=1:nargin(f)
        err_out = err_out + nderiv(f,i,y).^2.*err(:,i).^2;
    end
    err_out = sqrt(err_out);
end

function out = nderiv(f,i,x,h)
if nargin<4
    h = eps*1e6;
end
    x1 = mat2cell(x,size(x,1),ones(1,size(x,2)));
    h2 = zeros(size(x));
    h2(:,i) = h;
    x2 = mat2cell(x+h2,size(x,1),ones(1,size(x,2)));
    out = (f(x2{:})-f(x1{:}))/h;
end