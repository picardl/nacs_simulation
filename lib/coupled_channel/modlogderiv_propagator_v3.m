
function [Xout,Mout,node] = modlogderiv_propagator_v3(Xin,Wabc,h,E,m,back,diag_boo)
if nargin<7
    diag_boo = 1;
end

% array sizes
N = size(Xin,1);
I = eye(N);
notI = 1-I;

% Improve accuracy by changing basis so Wc is mostly diagonal
% Does this make sense? Seems to give very different results
% "Coupled channel methods for solving the bound-state Schrodinger 
% equation" by Jeremy Hutson says yes, this helps and is a known method.
if diag_boo
    [V,~] = eig(Wabc(:,:,2));
    Wabc(:,:,1) = V'*(Wabc(:,:,1)*V);
    Wabc(:,:,2) = V'*(Wabc(:,:,2)*V);
    Wabc(:,:,3) = V'*(Wabc(:,:,3)*V);
    Xin = V'*Xin*V;
end

% potential matrix at grid points
W1a = Wabc(:,:,1).*notI;
Wc = Wabc(:,:,2);
W1b = Wabc(:,:,3).*notI;
if ~diag_boo
    W1c = Wc.*notI;
    hsq6 = 6/h^2;
    W1til = (I/(I-W1c/hsq6) - I)*hsq6;
end

% wavevector
ksq = 2*m*(diag(Wc)-E);
kpos = (ksq >= 0);
kneg = ~kpos;
kabs = sqrt(ksq.*(kpos-kneg));

% solution to zeroth order potential
kp = kabs(kpos);
km = kabs(kneg);
kph = kp*h;
kmh = km*h;
A = 0*ksq;
B = A;
A(kpos) = kp.*coth(kph);
A(kneg) = km.*cot(kmh);
B(kpos) = kp.*csch(kph);
B(kneg) = km.*csc(kmh);
A = diag(A);
B = diag(B);

% propagation matrices Y_i(a,c), Y_i(c,b), i=1..4
h3 = h/3;
Y1ac = A + h3*W1a;
Y23 = B;
if diag_boo
    Y4ac_1cb = A;
else
    Y4ac_1cb = A + 2*h3*W1til; 
end
Y4cb = A + h3*W1b;

% full step propagation matrices
Zacb_inv = 2*Y4ac_1cb;
Y23ab = Y23*(Zacb_inv\Y23);
Y1ab = Y1ac - Y23ab;
Y4ab = Y4cb - Y23ab;

% propagate either backward or forward
if ~back
    Xa = Xin; % log-derivative wavefunction at point a
    Zab = Xa+Y1ab;
    Mout = Zab\Y23ab;
    Xb = Y4ab - Y23ab*Mout;
    Xout = Xb; % log-derivative wavefunction at point b
    
    Zab_eig = eig(Zab);
    node = sum(Zab_eig<0);
else
    Xb = Xin; % log-derivative wavefunction at point b
    Zab = Xb-Y4ab;
    Mout = -Zab\Y23ab;
    Xa = -Y1ab + Y23ab*Mout;
    Xout = Xa; % log-derivative wavefunction at point a
    
    Zab_eig = eig(Zab);
    node = sum(Zab_eig>0);
end

% return to original basis
if diag_boo
    Mout = (V*Mout)*V';
    Xout = (V*Xout)*V';
end

end