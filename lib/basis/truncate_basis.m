function [basis,rc_keep,Nchn] = truncate_basis(basis,op_fun,values,tol)


if nargin<4
    tol = 1e-3;
end

f = fields(basis);

if all(ismember(f,{'qnums','ops'}))
    no_sub_bases_boo = 1;
    b = basis;
    basis = struct();
    basis.temp = b;
    f = fields(basis);
else
    no_sub_bases_boo = 0;
end

f(strcmp(f,'change')) = [];

for i = 1:numel(f)
    Nstates = size(basis.(f{i}).qnums,1);
    [row,col] = ndgrid(1:Nstates,1:Nstates);
    
    op_test = op_fun(basis.(f{i}).ops);
    op_test_diag = diag(op_test);
    op_test_offdiag = op_test - diag(op_test_diag);
    
    i_keep = [];
    for j = 1:numel(values)
        i_keep = cat(1,i_keep,find(abs(diag(op_test) - values(j)) < tol));
    end
    i_keep = sort(i_keep);
    
    nz_offdiag = find(abs(triu(op_test_offdiag))>tol);
    [r_od,c_od] = ind2sub([Nstates Nstates],nz_offdiag);
    
    if mean(abs(diag(op_test_offdiag(r_od(~ismember(r_od,i_keep)),c_od(~ismember(r_od,i_keep))))))>tol
        warning(['basis ' f{i} ' could not be truncated cleanly. skipping.']);
        continue
    end
    
    Nchn.(f{i}) = numel(i_keep);
    rc_keep.(f{i}) = find(ismember(row,i_keep) & ismember(col,i_keep));
    
    g = fields(basis.(f{i}).ops);
    for j = 1:numel(g)
        basis.(f{i}).ops.(g{j}) = reshape(basis.(f{i}).ops.(g{j})(rc_keep.(f{i})),Nchn.(f{i}),Nchn.(f{i}));
    end
    
    basis.(f{i}).qnums = basis.(f{i}).qnums(i_keep,:);
    
    if ismember('change',fields(basis))
        c = fields(basis.change);
        fi_first = find(cellfun(@any,regexp(c,['^' f{i}],'match','once')))';
        fi_last = find(cellfun(@any,regexp(c,[f{i} '$'],'match','once')))';
        for j = fi_first
            basis.change.(c{j}) = basis.change.(c{j})(i_keep,:);
        end
        for j = fi_last
            
            basis.change.(c{j}) = basis.change.(c{j})(:,i_keep);
        end
        
    end
end

if no_sub_bases_boo
    basis = basis.temp;
end

end