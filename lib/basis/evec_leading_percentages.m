function ind = evec_leading_percentages(V,N)

pct = abs(V).^2;

ind = zeros(N,size(V,2));
for i = 1:N
    [~,ind(i,:)] = max(pct,[],1);
    pct(ind(i,:),:) = -inf;
end

end