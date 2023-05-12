function out = spin_rotation_element(r,c)

out = (-1).^(r.J+r.S-2*r.Omega) .* ...
    sqrt(r.J.*(r.J+1).*(2*r.J+1).*r.S.*(r.S+1).*(2*r.S+1)) .*...
    (w3j(r.J,-r.Omega,1,1,c.J,c.Omega).*w3j(r.S,-r.Omega,1,1,c.S,c.Omega) + ...
    w3j(r.J,-r.Omega,1,0,c.J,c.Omega).*w3j(r.S,-r.Omega,1,0,c.S,c.Omega) + ...
    w3j(r.J,-r.Omega,1,-1,c.J,c.Omega).*w3j(r.S,-r.Omega,1,-1,c.S,c.Omega))...
    .*(r.J==c.J).*(r.S==c.S).*(r.m_J==c.m_J);

end