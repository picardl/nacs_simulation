function ops = build_operators(basis)

ops = struct();

qnum_names = basis.Properties.VariableNames;

% identify j,m pairs
[a,b] = ismember(strcat('m_',qnum_names),qnum_names);
j_names = qnum_names(a);
m_names = qnum_names(b(a));

for i = 1:numel(j_names)
    J = j_names{i};
    mJ = m_names{i};
    ops.([J '_p']) = operator_matrix(@spher_op,basis,{J,mJ},1,J,mJ);
    ops.([J '_m']) = operator_matrix(@spher_op,basis,{J,mJ},-1,J,mJ);
    ops.([J '_x']) = (ops.([J '_m'])-ops.([J '_p']))/sqrt(2);
    ops.([J '_y']) = -(ops.([J '_m'])+ops.([J '_p']))/(sqrt(2)*1i);
    ops.([J '_z']) = operator_matrix(@spher_op,basis,{J,mJ},0,J,mJ);
    ops.([J '_sq']) = ops.([J '_x'])^2 + ops.([J '_y'])^2 + ops.([J '_z'])^2;
end
ops.I = eye(size(basis,1));


end