
function P = jacobi(n,alpha,beta,z)

prefactor = gamma(alpha+n+1)./(factorial(n).*gamma(alpha+beta+n+1));

x = reshape(0:n,1,1,n+1);

P = prefactor.*sum(cell2mat(arrayfun(@(m) nchoosek(n,m).*gamma(alpha+beta+n+m+1).*((z-1)/2).^(m)./gamma(alpha+m+1),x,'un',0)),3);

end