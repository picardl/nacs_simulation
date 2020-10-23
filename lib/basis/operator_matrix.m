function out = operator_matrix(element_fun,qnums,include,varargin)

if isa(qnums,'cell')
    qnums_row = qnums{1};
    qnums_col = qnums{2};
elseif isa(qnums,'table')
    qnums_row = qnums;
    qnums_col = qnums;
else
    error('qnums must be a cell array or a table');
end

Nr = size(qnums_row,1);
Nc = size(qnums_col,1);

[r,c] = ndgrid(1:Nr,1:Nc);

exclude_row = qnums_row.Properties.VariableNames(~ismember(qnums_row.Properties.VariableNames,include));
exclude_col = qnums_col.Properties.VariableNames(~ismember(qnums_col.Properties.VariableNames,include));

[row_in_col,row_col_ind] = ismember(exclude_row,exclude_col);
row_col_ind(row_col_ind==0) = [];
if any(~row_in_col)
    warning('columns have spectator quantum numbers that do not appear in rows')
end

boo = all(qnums_row{r,exclude_row} == qnums_col{c,exclude_col(:,row_col_ind)},2);

incl_row = include(ismember(include,qnums_row.Properties.VariableNames));
incl_col = include(ismember(include,qnums_col.Properties.VariableNames));
out = reshape(boo.*element_fun(qnums_row(r,incl_row),qnums_col(c,incl_col),varargin{:}),Nr,Nc,[]);

end