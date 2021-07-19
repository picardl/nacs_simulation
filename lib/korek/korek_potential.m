function out = korek_potential(r,term,index)

if nargin<1
    r = linspace(4.5,50,100);
end
if nargin<2
    term = '1S';
end
term = cellstr(term);
if nargin<3
    index = num2cell(ones(1,numel(term)));
elseif isnumeric(index)
    index = num2cell(index);
end

out = zeros(numel(r),numel(term));
for i = 1:numel(term)
    if isnan(str2double(term{i}(1)))
        index{i} = term{i}(1);
        term{i} = term{i}(2:end);
    end
    
    if ischar(index{i})
        if strcmp(index{i},'X')
            index{i} = 1;
        else
            alph = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
            singlet_ind = strfind(upper(alph),index{i});
            triplet_ind = strfind(lower(alph),index{i});
            if any(singlet_ind) && str2double(term{i}(1))==1
                index{i} = singlet_ind+1;
            elseif any(triplet_ind) && str2double(term{i}(1))==3
                index{i} = triplet_ind;
            else
                error('term symbol does not exist');
            end
        end
    end
    
    data = readtable(['NaCs' term{i} '.csv'],'headerlines',1);
    out(:,i) = interp1(data{:,1},data{:,index{i}+1},r(:),'spline','extrap');
end

end