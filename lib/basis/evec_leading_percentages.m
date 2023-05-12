function ind = evec_leading_percentages(V,N)

% Returns the row indices of the N leading percentages of the columns of
% matrix V. Each column of V is supposed to be a state vector psi_i of a
% quantum spin system, represented in some basis.

pct = abs(V).^2;

ind = zeros(N,size(V,2));
for i = 1:N
    [~,ind(i,:)] = max(pct,[],1);
    pct(ind(i,:),:) = -inf;
end

end