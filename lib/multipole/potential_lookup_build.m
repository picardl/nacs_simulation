function potential_lookup_build

n = 0:8;

powers = [];
coeffs = [];
l = [];
m = [];

for i = 1:numel(n)
    data = load(['mp' num2str(n(i)) '.txt']);
    L = data(1,4:end);
    M = data(2,4:end);
    powers = [powers; data(3:end,1:3)];
    coeffs{i+1} = data(3:end,4:end);
    l = [l, L];
    m = [m, M];
end

lm = [l;m]';
coeffs = blkdiag(coeffs{:});

save('potential_lookup','lm','powers','coeffs')

end