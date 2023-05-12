function basis = balance_operators(basis,b1,b2)


f = setdiff(fields(basis.(b1).ops),fields(basis.(b2).ops));
for i = 1:numel(f)
    basis.(b2).ops.(f{i}) = pagemtimes(basis.change.([b1 '_' b2]),'none', ...
        pagemtimes( basis.(b1).ops.(f{i}),'none' , basis.change.([b1 '_' b2]),'ctranspose'),'none');
end

f = setdiff(fields(basis.(b2).ops),fields(basis.(b1).ops));
for i = 1:numel(f)
    basis.(b1).ops.(f{i}) = pagemtimes(basis.change.([b1 '_' b2]),'ctranspose', ...
        pagemtimes( basis.(b2).ops.(f{i}),'none' , basis.change.([b1 '_' b2]),'none'),'none');
end

end