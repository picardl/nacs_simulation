function D = fourier_deriv_matrix(x,order)
% function D = FOURIER_DERIV_MATRIX(x,order) computes matrix D for taking 
% derivatives of vector x by Fourier transform methods. The input x must be 
% either a row or a column vector. 
% WBC 200323

N = numel(x);
k = fourier_fvec(x,1);

sz = size(x);
if sz(1)==1 && sz(2)>1
    % row vector
    dim = 2;
    sz_flag = 0;
elseif sz(1)>1 && sz(2)==1
    % column vector
    dim = 1;
    sz_flag = 0;
else
    sz_flag = 1;
end

if numel(sz)>2 || sz_flag
    error('x must be a row or column vector');
end

D = ifft((k.^order).*fft(eye(N),[],dim),[],dim);
D = (-1).^order * (1i.^order) * (D + D')/2;

end