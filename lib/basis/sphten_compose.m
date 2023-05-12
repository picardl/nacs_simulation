function out = sphten_compose(k,T1,T2)

% compose two spherical tensors of rank k1,p1 and k2,p2 into a tensor of
% rank k,p

k1 = (numel(T1)-1)/2;
k2 = (numel(T2)-1)/2;

if ~((abs(k1-k2)<=k) && (k<=(k1+k2)))
    error('cannot compose a tensor of rank %d from tensors of rank %d and rank %d',k,k1,k2);
end

if isa(T1,'sym') || isa(T2,'sym')
    out = sym(zeros(1,2*k+1));
else
    out = zeros(1,2*k+1);
end

% clebsch_table = array2table(zeros(0,7),'variablenames',{'k1','p1','k2','p2','k','p','CG'});
% ind = 1;
for p = -k:k
    for p1 = -k1:k1
        i1 = p1+k1+1;
        for p2 = -k2:k2
            i2 = p2+k2+1;
            out(p+k+1) = out(p+k+1) + clebsch(k1,p1,k2,p2,k,p) * T1(i1) * T2(i2);
%             clebsch_table{ind,:} = [k1,p1,k2,p2,k,p,clebsch(k1,p1,k2,p2,k,p)];
%             ind = ind+1;
        end
    end
end

end