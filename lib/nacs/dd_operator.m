function Vdd = dd_operator(ops,s1,s2)

s1tot = (ops.([s1 '_x'])+ops.([s1 '_y'])+ops.([s1 '_z']));
s2tot = (ops.([s2 '_x'])+ops.([s2 '_y'])+ops.([s2 '_z']));
Vdd = op_dot(ops,s1,s2) - 3*s1tot.*s2tot;

end