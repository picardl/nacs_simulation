
function w = w6j(varargin)
% FUNCTION W6J computes Wigner 6j coefficients. Mostly from the file
% exchange, modified by WBC to take matrix arguments.

if nargin==1
    x = varargin{1};
elseif nargin==6
    sz = reshape(cell2mat(cellfun(@(x) size(x),varargin,'un',0)),2,[])';
    varargin = arrayfun(@(i) repmat(varargin{i},max(sz)./sz(i,:)),1:nargin,'un',0);
    varargin = arrayfun(@(i) varargin{i}(:),1:nargin,'un',0);
    x = cell2mat(varargin);
end

w = arrayfun(@(i) Wigner6j(x(i,:)),1:size(x,1));
if nargin==6
    w = reshape(w,max(sz));
end

w = real(w);

    function w = Wigner6j(in)
        in = round(in*2)/2;
        
        % Finding Triangular coefficients
        
        tri1 = triangle_coeff(in(1),in(2),in(3));
        tri2 = triangle_coeff(in(1),in(5),in(6));
        tri3 = triangle_coeff(in(4),in(2),in(6));
        tri4 = triangle_coeff(in(4),in(5),in(3));
        
        if (tri1==0||tri2==0||tri3==0||tri4==0)
            w=0;
            return
        end
        
        % Finding the range of summation in the Racah formula.
        
        a(1) = in(1) + in(2) + in(3);
        a(2) = in(1) + in(5) + in(6);
        a(3) = in(4) + in(2) + in(6);
        a(4) = in(4) + in(5) + in(3);
        
        rangei = max(a);
        
        k(1) = in(1) + in(2) + in(4) + in(5);
        k(2) = in(2) + in(3) + in(5) + in(6);
        k(3) = in(3) + in(1) + in(6) + in(4);
        
        rangef = min(k);
        
        t = rangei:rangef;
        t = t(:)';
        
        w = sum(((-1).^t).*exp(gammaln(t+2) - fung(t,in(1),in(2),in(3),in(4),in(5),in(6))),2);
        
        w = sqrt(tri1*tri2*tri3*tri4)*w;
    end

    function r = fung(t, j1,j2,j3,J1,J2,J3)
        % Calculating the logarithm of the denominator in Racah Formula, using
        % the gamma function in place of the factorial.
        r = sum(gammaln([(t-j1-j2-j3);(t-j1-J2-J3);(t-J1-j2-J3);(t-J1-J2-j3);...
            (j1+j2+J1+J2-t);(j2+j3+J2+J3-t);(j3+j1+J3+J1-t)] + 1),1);
    end

    function tri = triangle_coeff(a,b,c)
        % Calculates triangle coefficients for angular momenta.
        % This version returns 0 if the triangle inequalities are violated.  (RAH)
        
        if (a<0 || b<0 || c<0)
            tri = 0;
            return
        end
        
        for xa = abs(a-b):(a+b)
            if c==xa
                tri = exp(sum(bsxfun(@times,[1;1;1;-1],gammaln([a+b-c;a-b+c;-a+b+c;a+b+c+1]+1)),1));
                return
            end
        end
        tri = 0;
        
    end

end