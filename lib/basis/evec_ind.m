function [qnums_ind,V_ind] = evec_ind(qnum_names,qnum_vals,basis,V)

if nargin<4
    V = eye(size(basis.qnums,1));
end

% get indices of a state with particular quantum numbers in a basis, and
% the column whose dominant admixture is that state in a matrix of
% eigenvectors

qnums_ind = find(all(abs(basis.qnums{:,qnum_names}-qnum_vals)<1e-10,2));

[~,V_ind] = max(abs(V(:,:,end)).^2,[],1);

V_ind = find(ismember(V_ind,qnums_ind));

end