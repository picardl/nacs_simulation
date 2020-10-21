function f = fourier_fvec(t,shift_boo)
% function f = FOURIER_FVEC(t,shift_boo)
% Fourier space f-vector for a given real space t-vector. The corresponding
% transformed Y-vector is obtained by fft(y). The input t must be either a
% row or a column vector. If shift_boo==1, the frequency vector is returned
% in the standard Matlab frequency ordering, while if shift_boo==0, the
% frequency vector is returned from minimum to maximum. 
% WBC 20200323

if nargin<2
    shift_boo = 1;
end

sz = size(t);

if sz(1)==1 && sz(2)>1
    % row vector
    t = t(:);
    col_boo = 0;
    sz_flag = 0;
elseif sz(1)>1 && sz(2)==1
    % column vector
    col_boo = 1;
    sz_flag = 0;
else
    sz_flag = 1;
end

if numel(sz)>2 || sz_flag
    error('x must be a row or column vector');
end

N = numel(t);
a = t(1);
b = t(end) + peak2peak(t)/(N-1);

dt_avg = mean(diff(t));
dt = t(2)-t(1);

if (dt-dt_avg)/dt_avg > 1e-6
    error('x must be equally spaced');
end

if mod(N,2)==0
    f = 2*pi/(b-a)*((-N/2):(N/2-1));
else
    f = 2*pi/(b-a)*(-floor(N/2):ceil(N/2-1));
end

if shift_boo
    f = ifftshift(f);
end

if col_boo
    f = f';
end

end