

a.qnums = build_basis({'J','S','Lambda'},{0:1,1,0},[1 1 1],'a');
a.ops = build_operators(a.qnums);

b.qnums = build_basis({'N','S','Lambda'},{0:2,1,0},[1 1 1],'b');
b.ops = build_operators(b.qnums);

b_a = operator_matrix(@case_b2a_element,{a.qnums,b.qnums},{'J','Omega','S','Sigma','N','Lambda'});
f = fields(b.ops);
for i = 1:numel(f)
    a.ops2.(f{i}) = b_a*b.ops.(f{i})*b_a';
    
    err(i) = max(reshape(abs(a.ops2.(f{i}) - a.ops.(f{i})),[],1));
end

err

f = fields(a.ops);
for i = 1:numel(f)
    b.ops2.(f{i}) = b_a'*a.ops.(f{i})*b_a;
    
    err(i) = max(reshape(abs(b.ops2.(f{i}) - b.ops.(f{i})),[],1));
end

err
