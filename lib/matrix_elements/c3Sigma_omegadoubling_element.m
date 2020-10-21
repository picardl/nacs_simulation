function out = c3Sigma_omegadoubling_element(r,c)

out = (r.Omega == -c.Omega).*(r.J==c.J).*r.J.*(r.J+1)/4;

end