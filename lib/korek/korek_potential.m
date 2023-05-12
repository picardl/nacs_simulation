function out = korek_potential(r,term)

if nargin<1
    r = linspace(4.5,50,100);
end

term = cellstr(term);

data = readtable('korek_data.xlsx');

out = interp1(data.R,data{:,term},r(:),'linear','extrap');

% out = zeros(numel(r),numel(term));
% for i = 1:numel(term)
%     pp = csape(data.R,data.(term{i}),'variational');
%     out(:,i) = fnval(pp,r);
% end

end