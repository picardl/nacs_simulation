function cj = clebsch(j1,m1,j2,m2,J,M)

cj = (-1).^(-j1+j2-M).*w3j(j1,m1,j2,m2,J,-M).*sqrt(2*J+1);

end