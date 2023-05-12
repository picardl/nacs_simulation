function psi = randbloch(N)
% random two-level state vectors by point picking on bloch sphere.

if nargin<1
    N = 1;
end

cart = normrnd(0,1,3,N);
normcart = bsxfun(@rdivide,cart,sqrt(sum(cart.^2,1)));

phi = atan2(normcart(2,:),normcart(1,:));
theta = acos2(normcart(3,:),1);

psi = [cos(theta/2); exp(1i*phi).*sin(theta/2)];

end