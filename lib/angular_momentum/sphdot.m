function C = sphdot(A,B)
% spherical tensor dot product

if numel(A)==3 && numel(B)~=3 && size(B,3)==3
    flag = 1;
elseif numel(B)==3 && numel(A)~=3 && size(A,3)==3
    t = A;
    A = B;
    B = t;
elseif size(A,3)==3 && size(B,3)==3
    flag = 0;
    pref = reshape([-1 1 -1],1,1,3);
    A = flip(A,3);
    C = sum(pref.*pagemtimes(A,B),3);
end

if flag
    A = reshape(A,1,1,3);
    pref = reshape(-1:1,1,1,3);
    C = sum(pref.*(A.*B),3);
end

end