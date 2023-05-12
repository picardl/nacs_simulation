function field_lookup_build

data = load('multipole_potential_lookup');

lm = data.lm;

[powers,coeffs{1}] = polydernq(data.powers,data.coeffs,1);
[~,coeffs{2}] = polydernq(data.powers,data.coeffs,2);
[~,coeffs{3}] = polydernq(data.powers,data.coeffs,3);

coeffs = cellfun(@(x) -x,coeffs,'un',0);
coeffs = cat(3,coeffs{:});

save('field_lookup','powers','coeffs','lm')

end