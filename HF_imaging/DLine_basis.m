function basis = DLine_basis(atom,DLine)

c = constants();

name = [atom,DLine];

if nargin<1
    name = 'CsD2';
end

%Define quantum numbers and build basis and fundamental angular momentum
%operators
i_name = ['i_' name];
j_name = ['j_' name];
basis.qnums = build_basis({i_name,j_name},{c.(i_name),c.(j_name)},[1 1]);
basis.ops = build_operators(basis.qnums);

% Hyperfine Hamiltonian (see Steck alkali data)
eyeQ = eye(size(basis.qnums,1));
Ham_HF = c.(['zeta_hf_' name])*op_dot(basis.ops,i_name,j_name);
B_HF = c.(['B_hf_' name])*(3*op_dot(basis.ops,i_name,j_name)^2 + (3/2)*op_dot(basis.ops,i_name,j_name)...
    - eyeQ.*c.(j_name)*(c.(j_name)+1)*c.(i_name)*(c.(i_name) + 1))/(2*c.(i_name)*(2*c.(i_name)-1)*c.(j_name)*(2*c.(j_name)-1));
C_HF = c.(['C_hf_' name])*(10*op_dot(basis.ops,i_name,j_name)^3 + 20*op_dot(basis.ops,i_name,j_name)^2 + 2*op_dot(basis.ops,i_name,j_name)*...
    (eyeQ.*c.(i_name)*(c.(i_name)+1) + eyeQ.*c.(j_name)*(c.(j_name)+1)+ 3*eyeQ - eyeQ*3.*c.(i_name)*(c.(i_name)+1).*c.(j_name)*(c.(j_name)+1))...
    - eyeQ.*5.*c.(i_name)*(c.(i_name)+1).*c.(j_name)*(c.(j_name)+1));

TE = vec2sphTen([sin(theta)*cos(gamma),sin(theta)*1i*sin(gamma),cos(theta)],[sin(theta)*cos(gamma),-sin(theta)*1i*sin(gamma),cos(theta)]);
basis.UC.ops.LS = operator_matrix(@atom_ls,basis.UC.qnums,{'N','S','m_N'},const.X1Sigma.a_perp,const.X1Sigma.a_par,TE);

if DLine == 'D2' || DLine == 'D1' 
    Ham_HF = Ham_HF + B_HF + C_HF;
end
basis.ops.(['H_' name '_hyperfine']) = Ham_HF;

% Zeeman Hamiltonian
basis.ops.(['H0_' name '_zeeman']) = (c.(['gj_' name])*c.uB*basis.ops.([j_name '_z']) + c.(['gi_' atom])*c.uB*basis.ops.([i_name '_z']));

end
