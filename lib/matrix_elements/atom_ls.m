function out = atom_ls(r,c,TE)

pol_lab = zeros(3,5);
out = zeros(height(r),1);
for i = 1:length(out)
    for k = 1:3
        for p = 3-(k-1):3+(k-1)
            pol_lab(k,p) =  pol_lab(k,p)*wigner3j(r.J,-r.m_J,k,p,c.J,-c.m_J);
        end
    end
    out(i) = -(0.5)*spher_dot(pol_lab,TE);
end

end
