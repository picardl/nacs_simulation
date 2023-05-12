function q = basis_case_a(Lambda,S,Jrange,Omega_abs)

Lambda = unique([-Lambda Lambda]);
Sigma = -S:S;

[Lambda,Sigma] = ndgrid(Lambda,Sigma);
Lambda = Lambda(:);
Sigma = Sigma(:);
Omega = Lambda + Sigma;
S = S*ones(size(Omega));

J = 0:max(Jrange);
m_J = -max(J):max(J);

[J,m_J,p] = ndgrid(J,m_J,1:numel(Omega));
J = J(:);
m_J = m_J(:);
p = p(:);

q = table();

q.S = S(p);
q.Lambda = Lambda(p);
q.Sigma = Sigma(p);
q.Omega = Omega(p);
q.J = J;
q.m_J = m_J;

del = abs(q.Omega)>q.J | q.J<min(Jrange) | abs(q.m_J)>q.J;
q(del,:) = [];

if nargin==4
    del = abs(q.Omega)~=Omega_abs;
    q(del,:) = [];
end

end