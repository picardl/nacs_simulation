function out = case_b2a_element(r,c)

% rows are hund's case (a)

out = (-1).^(c.N-r.J+r.Sigma) .* clebsch(r.J,r.Omega,r.S,-r.Sigma,c.N,c.Lambda) .* (r.J==c.J).*(r.S==c.S);

end