function out = Sz_UC_element(r,c)

out = zeros(size(r.S));
for p = -1:1
    out = out + (-1).^p .* (-1).^(r.S-r.m_S) .* w3j(r.S,-r.m_S,1,-p,c.S,c.m_S) .* sqrt(r.S.*(r.S+1).*(2*r.S+1)) .*(r.S==c.S) .*...
        (-1).^(r.m_N-r.Lambda) .* sqrt((2*r.N+1).*(2*c.N+1)).*w3j(r.N,-r.m_N,1,p,c.N,c.m_N).*w3j(r.N,-r.Lambda,1,0,c.N,c.Lambda);
end

end