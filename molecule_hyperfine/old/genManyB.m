Bs = [850:1:860]*1e-4;
save_basis = 'aUC';
for b = 1:length(Bs)
    B = Bs(b);
    feshbach(B,save_basis,2)
%     c3Sigma(B,save_basis,2) 
end