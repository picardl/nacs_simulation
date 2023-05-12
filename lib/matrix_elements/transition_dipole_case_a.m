function out = transition_dipole_case_a(r,c,p,T1q)

if nargin==4
    out = (-1).^(r.m_J-r.Omega) .* sqrt((2*r.J+1).*(2*c.J+1)) ...
        .*w3j(r.J,-r.m_J,1,p,c.J,c.m_J).*...
        (w3j(r.J,-r.Omega,1,1,c.J,c.Omega).*T1q(1) ...
        + w3j(r.J,-r.Omega,1,0,c.J,c.Omega).*T1q(2) ...
        + w3j(r.J,-r.Omega,1,-1,c.J,c.Omega).*T1q(3));
else
    out = (-1).^(r.m_J-r.Omega) .* sqrt((2*r.J+1).*(2*c.J+1)) ...
        .*w3j(r.J,-r.m_J,1,p,c.J,c.m_J).*w3j(r.J,-r.Omega,1,r.Omega-c.Omega,c.J,c.Omega);
end

if ismember('eta',r.Properties.VariableNames)
    boo = abs(r.eta-c.eta)==1;
    out = out.*boo;
end

end