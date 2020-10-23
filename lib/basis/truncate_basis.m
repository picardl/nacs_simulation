function [basis,rc_keep,Nchn] = truncate_basis(basis,op_fun,values,tol)

if nargin<4
    tol = 1e-3;
end

f = fields(basis);

for i = 1:numel(f)
    Nstates = size(basis.(f{i}).qnums,1);
    [row,col] = ndgrid(1:Nstates,1:Nstates);
    
    i_keep = [];
    for j = 1:numel(values)
        i_keep = cat(1,i_keep,find(abs(diag(op_fun(basis.(f{i}).ops)) - values(j)) < tol));
    end
    
    
    Nchn.(f{i}) = numel(i_keep);
    rc_keep.(f{i}) = find(ismember(row,i_keep) & ismember(col,i_keep));
    
    g = fields(basis.(f{i}).ops);
    for j = 1:numel(g)
        basis.(f{i}).ops.(g{j}) = reshape(basis.(f{i}).ops.(g{j})(rc_keep.(f{i})),Nchn.(f{i}),Nchn.(f{i}));
    end
    
    basis.(f{i}).qnums = basis.(f{i}).qnums(i_keep,:);
    
end

end