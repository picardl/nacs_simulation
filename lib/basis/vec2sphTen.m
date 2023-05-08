function T = vec2sphTen(v1,v2)
%VEC2SPHTEN Convert 3 element vector to spherical tensor indexed as T(k,p)
T = zeros(3,5);
if nargin == 1
    T(2,3) = vec(3);
    T(2,2) = (1/sqrt(2))*(vec(1) - 1i*vec(2));
    T(2,4) = -(1/sqrt(2))*(vec(1) + 1i*vec(2));
else
    cr = cross(v1,v2);
    elw = v1.*v2;
    T(1,3) = (-1/sqrt(3))*sum(elw);
    
    T(2,3) = 1i/sqrt(2)*cr(3);
    T(2,2) = 1/sqrt(2)*(-1i*cr(1) - cr(2));
    T(2,4) = 1/sqrt(2)*(1i*cr(1) - cr(2));
    
    T(3,3) = 1/sqrt(6)*(2*elw(3) - elw(1) - elw(2));
    T(3,2) = 1/2*((v1(1)*v2(3) + v1(3)*v2(1)) + 1i*(v1(2)*v2(3) + v1(3)*v2(2)));
    T(3,4) = 1/2*(-(v1(1)*v2(3) + v1(3)*v2(1)) + 1i*(v1(2)*v2(3) + v1(3)*v2(2)));
    T(3,1) = 1/2*(elw(1) - elw(2) - 1i*(v1(1)*v2(2)+v1(2)*v2(1)));
    T(3,5) = 1/2*(elw(1) - elw(2) + 1i*(v1(1)*v2(2)+v1(2)*v2(1)));
end


