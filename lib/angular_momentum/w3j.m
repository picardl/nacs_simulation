
function w = w3j(varargin)
% FUNCTION W3J computes Wigner 3j coefficients. Mostly from the file
% exchange, modified by WBC to take matrix arguments.

%% process arguments
if nargin==1
    x = varargin{1};
elseif nargin==6
    sz = reshape(cell2mat(cellfun(@(x) size(x),varargin,'un',0)),[],nargin)';
    varargin = arrayfun(@(i) repmat(varargin{i},max(sz)./sz(i,:)),1:6,'un',0);
    varargin = arrayfun(@(i) varargin{i}(:),1:6,'un',0);
    x = cell2mat(varargin);
end

x = round(x*2)/2;

%% calculate only unique elements
[x,~,reverse_unique_ind] = unique(x,'rows');

%% use selection rules to find zero 3j's
j1 = x(:,1);
m1 = x(:,2);
j2 = x(:,3);
m2 = x(:,4);
j3 = x(:,5);
m3 = x(:,6);

j123 = [j1 j2 j3];
m123 = [m1 m2 m3];

unphys = any(j123<0 + rem(j123-m123,1),2) + any(rem(x,1/2),2);
selection = (j3>(j1+j2)) + (j3<abs(j1-j2)) + (abs(sum(m123,2))>0) + any(abs(m123)>j123,2);
simple = bsxfun(@and,~any(m123,2),rem(sum(j123,2),2));

zero = bsxfun(@or,bsxfun(@or,unphys,selection),simple);

%% calculate nonzero 3j's
x(zero,:) = [];

j1 = x(:,1);
m1 = x(:,2);
j2 = x(:,3);
m2 = x(:,4);
j3 = x(:,5);
m3 = x(:,6);

t1 = j2 - m1 - j3;
t2 = j1 + m2 - j3;
t3 = j1 + j2 - j3;
t4 = j1 - m1;
t5 = j2 + m2;

T = arrayfun(@(a,b,c,d,e) max(0,max(a,b)):min(c,min(d,e)),t1,t2,t3,t4,t5,'un',0);
temp = arrayfun(@(i) f(j1(i),j2(i),j3(i),m1(i),m2(i),m3(i),T{i}),1:length(T));

%% assign and reshape output
w(~zero) = temp;
w(zero) = 0;
w = w(reverse_unique_ind);

if nargin==6
    w = reshape(w,max(sz));
end

w = real(w);

    function y = f(j1,j2,j3,m1,m2,m3,t)
        t1 = j2 - m1 - j3;
        t2 = j1 + m2 - j3;
        t3 = j1 + j2 - j3;
        t4 = j1 - m1;
        t5 = j2 + m2;
        y = sum( (-1).^t .* exp( -ones(1,6) * gammaln( [t; t-t1; t-t2; t3-t; t4-t; t5-t] +1 ) + ...
            gammaln( [j1+j2+j3+1, j1+j2-j3, j1-j2+j3, -j1+j2+j3, j1+m1, j1-m1, j2+m2, j2-m2, j3+m3, j3-m3] +1 ) ...
            * [-1; ones(9,1)] * 1/2 ) ) * (-1)^( j1-j2-m3 );
    end

end