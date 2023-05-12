function x = herm(x)

x = (x + conj(permute(x,[2 1 3:numel(size(x))])))/2;

end