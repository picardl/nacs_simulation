function out = mat_el_nd(A,i,j)

sz = size(A);
if numel(sz)<3
    out = A(i,j);
else
    out = A(sub2ind(sz,i*ones(size(1:sz(3))),j*ones(size(1:sz(3))),1:sz(3)));
end


end