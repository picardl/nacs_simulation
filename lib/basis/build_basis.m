function q = build_basis(fixed_qnum_names,fixed_qnum_values,m_boo,hunds_case)

grid_input = {};
m_cols = [];
k_cols = [];
all_qnum_names = {};

for i = 1:numel(fixed_qnum_names)
    j_name = fixed_qnum_names{i};
    
    ji = fixed_qnum_values{i};
    if m_boo(i)==1
        k_name = [];
        m_name = ['m_' j_name];
        mi = -max(ji):max(ji);
        ki = [];
        m_cols_boo = [0 1];
        k_cols_boo = [0 0];
    elseif m_boo(i)==2
        switch j_name
            case 'J'
                k_name = 'Omega';
            case 'N'
                k_name = 'Lambda';
            case 'S'
                k_name = 'Sigma';
            otherwise
                k_name = ['k_' j_name];
        end
        m_name = ['m_' j_name];
        ki = unique([-abs(min(ji)) abs(min(ji))]);
        mi = -max(ji):max(ji);
        m_cols_boo = [0 0 1];
        k_cols_boo = [0 1 0];
    else
        k_name = [];
        m_name = [];
        mi = [];
        ki = [];
        m_cols_boo = [0];
        k_cols_boo = [0];
    end
    
    grid_input = cat(2,grid_input,ji,ki,mi);
    m_cols = cat(2,m_cols,m_cols_boo);
    k_cols = cat(2,k_cols,k_cols_boo);
    all_qnum_names = cat(2,all_qnum_names,j_name,k_name,m_name);
end
num_qnums = numel(all_qnum_names);

[q{1:num_qnums}] = ndgrid(grid_input{:});
q = cell2mat(cellfun(@(x) x(:),q,'un',0));

q = sortrows(q);
q = array2table(q);
q.Properties.VariableNames = all_qnum_names;

m_names = all_qnum_names(m_cols>0);
del = [];
for i = 1:numel(m_names)
    del = cat(1,del,find(abs(q.(m_names{i})) > q.(m_names{i}(3:end))));
end
q(del,:) = [];

if nargin==4
    switch hunds_case
        case 'a'
            q.Sigma = q.Omega-q.Lambda;
        otherwise
            warning('You asked for hund''s case %s but it is not implemented.',hunds_case)
    end
end


end