function out = spher_dot(T1,T2)
%SPHER_DOT Calculate the dot product of two spherical tensor operators T1 and
%T2, with dimensions k, p.
%See https://bartmcguyer.com/notes/note-9-LightShifts.pdf
%% Check input sizes
k1 = size(T1,1);
k2 = size(T2,1);
if k1 == 2 && k2 == 3
    T2 = T2(1:2,:);
    disp('Truncating T2 to vector component only')
    kMax = 2;
elseif k1 == 3 && k2 == 2
    T1 = T1(1:2,:);
    disp('Truncating T1 to vector component only')
    kMax = 2;
elseif k1 ~= k2
    error('Input tensors must be 2x3 or 3x5')
else
    kMax = k1;
end
if size(T1,2) ~= 2*kMax - 1 || size(T2,2) ~= 2*kMax - 1
    error('Input tensors must be 2x3 or 3x5')
end

out = 0;
for k = 1:kMax
   for p = -(k-1):(k-1)
       out = out + (-1)^(k-1 - p)*0.5*T1(k,p+3)*T2(k,-p+3);
   end
end
%%

end

