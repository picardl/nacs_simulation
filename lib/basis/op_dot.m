
function j1dotj2 = op_dot(ops,j1_name,j2_name)

j1dotj2 = 0*ops.([j1_name '_x']);
for Q = {'_x','_y','_z'}
    q = [Q{:}];
    j1dotj2 = j1dotj2 + ops.([j1_name q])*ops.([j2_name q]);
end

end
