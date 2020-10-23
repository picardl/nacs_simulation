function basis = atom_basis(name)

c = constants();

if nargin<1
    name = 'Na';
end

i_name = ['i_' name];
s_name = ['s_' name];

basis.qnums = build_basis({i_name,s_name},{c.(i_name),c.(s_name)},[1 1]);
basis.ops = build_operators(basis.qnums);

% hyperfine
basis.ops.(['H_' name '_hyperfine']) = c.(['zeta_hf_' name])*op_dot(basis.ops,i_name,s_name);

% zeeman
basis.ops.(['H0_' name '_zeeman']) = (c.(['gs_' name])*c.uB*basis.ops.([s_name '_z']) + c.(['gi_' name])*c.uB*basis.ops.([i_name '_z']));

end