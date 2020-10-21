function b = join_basis(b1,b2)
% merge bases, qnums and operators using kronecker product

n1 = size(b1.qnums,1);
n2 = size(b2.qnums,1);

[i1,i2] = ndgrid(1:n1,1:n2);

i1 = i1(:);
i2 = i2(:);

b.qnums = [b1.qnums(i1,:) b2.qnums(i2,:)];

I1 = eye(n1);
I2 = eye(n2);

f1 = fields(b1.ops);
for i = 1:numel(f1)
    b.ops.(f1{i}) = kron(I2,b1.ops.(f1{i}));
end

f2 = fields(b2.ops);
for i = 1:numel(f2)
    b.ops.(f2{i}) = kron(b2.ops.(f2{i}),I1);
end

end