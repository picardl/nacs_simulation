function out = case_b2a_element(r,c)

% rows are hund's case (a)

% if you want to exactly translate from hund's case (b) to (a), you need to
% have one more N in case (b) than you have J in case(a).

% conversely to exactly translate from case (a) to case (b), you must have
% one more J in your case (a) basis than you have N in the case (b) basis

out = (-1).^(c.N-r.S+r.Omega) .*sqrt(2*c.N+1) .* w3j(r.J,r.Omega,r.S,-r.Sigma,c.N,-c.Lambda) .* (r.J==c.J).*(r.S==c.S);

end