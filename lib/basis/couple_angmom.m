
function [basis_cpl,basis_unc,basis_change_matrix] = couple_angmom(basis_unc,j1_name,j2_name,j3_name)

% couple angular momenta j1 and j2 to make new angular momentum j3=j1+j2,
% and change basis to one where j3^2 and j3z are diagonal.

qnums_unc = basis_unc.qnums;

basis_cpl = struct();
basis_cpl.qnums = couple_qnums(qnums_unc,j1_name,j2_name,j3_name);
basis_cpl.ops = struct();

m1_name = ['m_' j1_name];
m2_name = ['m_' j2_name];
m3_name = ['m_' j3_name];

% identify spectator quantum numbers
spec_qnum_names = qnums_unc.Properties.VariableNames(...
    ~ismember(qnums_unc.Properties.VariableNames,...
    {j1_name m1_name j2_name m2_name}));

Nstates = size(basis_unc.qnums,1);
[row,col] = ndgrid(1:Nstates,1:Nstates);
row = row(:);
col = col(:);

% create x,y,z,sq operators for coupled angular momentum in uncoupled basis
for Q = {'_x','_y','_z','_p','_m'}
    q = [Q{:}];
    basis_unc.ops.([j3_name q]) = basis_unc.ops.([j1_name q]) + basis_unc.ops.([j2_name q]);
end

basis_unc.ops.([j3_name '_sq']) = zeros(Nstates);
for Q = {'_x','_y','_z'}
    q = [Q{:}];
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
    basis_cpl.ops.(op_names{i}) = pagemtimes(basis_change_matrix,'ctranspose',pagemtimes(basis_unc.ops.(op_names{i}),basis_change_matrix),'none');
end

end