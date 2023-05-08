function out = mol_ls(r,c,a_perp,a_par,TE)

%molecule polarizability tensor
pol_mol = zeros(3,5);
pol_mol(1,3) = -(1/sqrt(3))*(2*a_perp + a_par);
pol_mol(3,3) = (1/sqrt(6))*(2*a_par -2*a_perp);

d_el = zeros(3,5,5,height(r));
for k = 1:3
    for p = 3-(k-1):3+(k-1)
        el = 0;
        for q = 3-(k-1):3+(k-1)
            d_el(k,p,q,:) = DMat(k-1,p-3,q-3,r.N,r.m_N,0,c.N,c.m_N,0);
        end
    end
end

pol_lab = zeros(3,5);
out = zeros(height(r),1);
for i = 1:length(out)
    for k = 1:3
        for p = 3-(k-1):3+(k-1)
            el = 0;
            for q = 1:5
                el = el + pol_mol(k,q)*d_el(k,p,q,i);
            end
            pol_lab(k,p) = el;
        end
    end
    out(i) = -(0.5)*spher_dot(pol_lab,TE);
end

out = out ./ (1/3 * (a_perp + 2*a_par));

function D_el = DMat(k,p,q,J1,m1,Omega1,J2,m2,Omega2)
 D_el = (-1).^(m1-Omega1).*sqrt((2*J1+1).*(2*J2+1)).*w3j(J1,-Omega1,k,q,J2,Omega2).*w3j(J1,-m1,k,p,J2,m2);
end

end