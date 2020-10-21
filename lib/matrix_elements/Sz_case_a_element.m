function out = Sz_case_a_element(r,c)

out = (-1).^(r.m_J-r.k_J) .* w3j(r.J,-r.m_J,1,0,c.J,c.m_J) .* w3j(r.J,-r.k_J,1,0,c.J,c.k_J) .* sqrt((2*r.J+1).*(2*c.J+1)) .* r.k_J;

end