
function [basis_cpl,basis_unc,basis_change_matrix] = couple_angmom(basis_unc,j1_name,j2_name,j3_name)

m1_name = ['m_' j1_name];
m2_name = ['m_' j2_name];
m3_name = ['m_' j3_name];

% identify spectator quantum numbers
spec_qnum_names = basis_unc.qnums.Properties.VariableNames(...
    ~ismember(basis_unc.qnums.Properties.VariableNames,...
    {j1_name m1_name j2_name m2_name}));

all_qnum_names = [j1_name j2_name j3_name m3_name spec_qnum_names];

% angular momenta to be coupled
uj1 = unique(basis_unc.qnums.(j1_name));
uj2 = unique(basis_unc.qnums.(j2_name));

% generate values of coupled quantum number
coupled_qnums = zeros(0,size(all_qnum_names,2));
for q = 1:numel(uj1)
    for r = 1:numel(uj2)
        spec_qnums_qr = unique(basis_unc.qnums{basis_unc.qnums.(j1_name)==uj1(q) & basis_unc.qnums.(j2_name)==uj2(r),spec_qnum_names},'rows');
        
        j3new = (abs(uj1(q)-uj2(r)):(uj1(q)+uj2(r)))';
        for s = 1:numel(j3new)
            m3 = (-j3new(s):j3new(s))'; 
            j3 = j3new(s)*ones(size(m3));
            j1 = uj1(q)*ones(size(m3));
            j2 = uj2(r)*ones(size(m3));
            
            coupled_qnums_qr = [j1 j2 j3 m3];
            
            [q_new_ind,q_spec_ind] = ndgrid(1:size(coupled_qnums_qr,1),1:size(spec_qnums_qr,1));
            q_new_ind = q_new_ind(:);
            q_spec_ind = q_spec_ind(:);
            
            coupled_qnums_qr = cat(2,coupled_qnums_qr(q_new_ind,:),spec_qnums_qr(q_spec_ind,:));
            
            coupled_qnums = cat(1,coupled_qnums,coupled_qnums_qr);
        end
    end
end

coupled_qnums = array2table(coupled_qnums,'VariableNames',all_qnum_names);

basis_cpl = struct();
basis_cpl.qnums = coupled_qnums;
basis_cpl.ops = struct();

Nstates = size(basis_unc.qnums,1);
[row,col] = ndgrid(1:Nstates,1:Nstates);
row = row(:);
col = col(:);

% create x,y,z,sq operators for coupled angular momentum in uncoupled basis
basis_unc.ops.([j3_name '_sq']) = zeros(Nstates);
for Q = {'_x','_y','_z'}
    q = [Q{:}];
    basis_unc.ops.([j3_name q]) = basis_unc.ops.([j1_name q]) + basis_unc.ops.([j2_name q]);
    basis_unc.ops.([j3_name '_sq']) = basis_unc.ops.([j3_name '_sq']) + basis_unc.ops.([j3_name q])^2;
end

% change basis to coupled one
basis_change_matrix = reshape(...
    clebsch(basis_unc.qnums{row,j1_name},basis_unc.qnums{row,m1_name},...
    basis_unc.qnums{row,j2_name},basis_unc.qnums{row,m2_name},...
    basis_cpl.qnums{col,j3_name},basis_cpl.qnums{col,m3_name}).*...
    all(basis_unc.qnums{row,[j1_name j2_name spec_qnum_names]}==basis_cpl.qnums{col,[j1_name j2_name spec_qnum_names]},2),...
    [Nstates Nstates]);

op_names = fields(basis_unc.ops);
for i = 1:numel(op_names)
    basis_cpl.ops.(op_names{i}) = basis_change_matrix'*(basis_unc.ops.(op_names{i})*basis_change_matrix);
end

end