function [a,a_b] = basis_b2a(b)

if any(~ismember({'N','S','J','Lambda'},b.qnums.Properties.VariableNames))
    error('basis provided is not Hunds case b')
end

b_qnums = unique(b.qnums(:,{'N','S','J','Lambda'}),'rows');

S = unique(b_qnums.S);
Lambda = unique(abs(b_qnums.Lambda));
Jmin = min(b_qnums.J);
Jmax = max(b_qnums.J);

b_qnums = unique(rmcol(b.qnums,{'N','J','m_J'}),'rows');

a_qnums = build_basis({'J','S','Lambda'},{Jmin:Jmax,S,Lambda},[0 0 0],'a');
a.qnums = outerjoin(b_qnums,a_qnums,'mergekeys',true);
a.ops = build_operators(a.qnums);
a_b = operator_matrix(@case_b2a_element,{a.qnums,b.qnums},{'J','Omega','S','Sigma','N','Lambda'});

% transform operators into case a
f = setdiff(fields(b.ops),fields(a.ops));
for i = 1:numel(f)
    a.ops.(f{i}) = pagemtimes(a_b,pagemtimes( b.ops.(f{i}),'none' , a_b,'ctranspose'));
end

end