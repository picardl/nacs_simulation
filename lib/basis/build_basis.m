function q = build_basis(fixed_qnum_names,fixed_qnum_values,m_boo,hunds_case)

% fixed_qnum_names: cell array of strings
% fixed_qnum_values: cell array of doubles or double arrays
% m_boo: double array with values 0,1, or 2.
%           0: don't create projection quantum numbers for this fixed qnum
%           1: do create projection quantum numbers (labeled m)
%           2: molecule case -- hund's case a,b,c. Create mol frame
%           projection quantum numbers with special names: Omega, Lambda,
%           Sigma.
% hunds_case: specify the hund's case for a molecular basis

%% molecule-specific basis
if nargin>3
    switch hunds_case
        case 'a'
            mol_qnum_names = {'J','Omega','S','Sigma','Lambda'};
        case 'b'
            mol_qnum_names = {'N','Lambda','S','J'};
        case 'c'
            mol_qnum_names = {'J','Omega'};
        otherwise
            error('that hunds case doesnt exist or hasnt been programmed');
    end
    
    [~,mol_qnum_ind] = ismember(mol_qnum_names,fixed_qnum_names);
    mol_qnum_ind = mol_qnum_ind(mol_qnum_ind~=0);
    mol_qnum_names = fixed_qnum_names(mol_qnum_ind(mol_qnum_ind~=0));
    mol_qnum_vals = fixed_qnum_values(mol_qnum_ind(mol_qnum_ind~=0));
    
    Lambda_Omega = mol_qnum_vals{...
        strcmp(mol_qnum_names,'Lambda') | (strcmp(mol_qnum_names,'Omega') & strcmp(hunds_case,'c'))};
    
    Nmax_Jmax = max(mol_qnum_vals{strcmp(mol_qnum_names,'N') | strcmp(mol_qnum_names,'J')});
    
    q_mol = mol_basis(hunds_case,Lambda_Omega,Nmax_Jmax,...
        mol_qnum_vals{strcmp(mol_qnum_names,'S')},...
        mol_qnum_vals{strcmp(mol_qnum_names,'Omega') & strcmp(hunds_case,'a')});
%     q_mol = mol_basis(hunds_case,Lambda_Omega,Nmax_Jmax,...
%         mol_qnum_vals{strcmp(mol_qnum_names,'S')});
    
    fixed_qnum_names(mol_qnum_ind) = [];
    fixed_qnum_values(mol_qnum_ind) = [];
    m_boo(mol_qnum_ind) = [];
end

%% vanilla basis
if ~isempty(fixed_qnum_names)
    grid_input = {};
    m_cols = [];
    all_qnum_names = {};
    
    
    for i = 1:numel(fixed_qnum_names)
        j_name = fixed_qnum_names{i};
        
        ji = fixed_qnum_values{i};
        if m_boo(i)==1
            m_name = ['m_' j_name];
            mi = -max(ji):max(ji);
            m_cols_boo = [0 1];
        else
            m_name = [];
            mi = [];
            m_cols_boo = [0];
        end
        
        grid_input = cat(2,grid_input,ji,mi);
        m_cols = cat(2,m_cols,m_cols_boo);
        all_qnum_names = cat(2,all_qnum_names,j_name,m_name);
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
    
    if nargin>3
        [i1,i2] = ndgrid(1:size(q,1),1:size(q_mol,1));
        q = [q(i1(:),:) q_mol(i2(:),:)];
    end
    
else
    q = q_mol;
end



end