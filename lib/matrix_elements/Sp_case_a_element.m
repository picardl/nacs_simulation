function out = Sp_case_a_element(r,c,p)

out = (-1).^(r.m_J-r.Omega) .* w3j(r.J,-r.m_J,1,p,c.J,c.m_J) .* w3j(r.J,-r.Omega,1,0,c.J,c.Omega) .* sqrt((2*r.J+1).*(2*c.J+1)) .* r.Omega;

end