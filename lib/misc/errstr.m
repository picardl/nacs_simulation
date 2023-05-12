
function str = errstr(m,s,d)
% ERRSTR    Prints a value and standard error in parenthesis.
%   ERRSTR(M,S) will print M(S), where the value of M is rounded
%   off to the precision of S.
%
%   Example: ERRSTR(12.77,0.13) will print 12.8(1)
%
%   ERRSTR(M,S,D) will print M(S) with D digits of precision.
%
%   Example: ERRSTR(12.77,0.13,2) will print 12.77(13)
%
%   ERRSTR(X) will print ERRSTR(MEAN(X),STDERR(X)).
if nargin<2
    s = stderr(m);
    m = mean(m);
end
if nargin<3
    d=1;
end

M=m(:);
S=s(:);
str = cell(1,length(M));
for i=1:length(M);
    m=M(i);
    s=S(i);
    
    if s==0
        str{i} = sprintf('%g',m);
    else
        x = 10^(floor(log10(s))+1-d);

        % deal with the corner case of s that should round up to 1. 
        if round(s)==1 && round(s/x)==10
            x = 1;
        end
        
        str{i} = sprintf(['%.' sprintf('%d',max(0,-log10(x))) 'f(%g)'],round(m/x)*x,round(s/x)*max(x,s<1));
    end
    
end % for

if length(M)==1
    str = str{i};
end

end % errstr