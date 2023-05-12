function out = diag_nd(A)

sz = size(A);
if sum(sz>1)==1
    out = A(:)';
elseif numel(sz)<3
    out = diag(A);
else
    if sz(1)~=sz(2)
        error('matrix must be size N*N*p');
    end
    out = A(bsxfun(@plus,(1:sz(1)+1:sz(1)^2)',(0:sz(3)-1)*sz(1)^2));
end

end