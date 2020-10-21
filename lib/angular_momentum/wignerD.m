
function D = wignerD(j,mp,m,alpha,beta,gamma)

k = min([j+m j-m j+mp j-mp]);

if k==j+m
    a = mp - m;
    lambda = mp - m;
elseif k==j-m
    a = m - mp;
    lambda = 0;
elseif k==j+mp
    a = m - mp;
    lambda = 0;
elseif k==j-mp
    a = mp - m;
    lambda = mp - m;
end

b = 2*j-2*k-a;

if (k+a > 2*j-k) || (b > k+b)
    d = 0*beta;
else
    d = (-1).^(lambda).*sqrt(nchoosek(2*j-k,k+a)./nchoosek(k+b,b)).*...
        sin(beta/2).^a.*cos(beta/2).^b.*jacobi(k,a,b,cos(beta));
end

D = exp(-1i*mp*alpha).*d.*exp(-1i*m*gamma);

end