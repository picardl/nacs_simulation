function X = rmcol(X,rmvars)
% funciton rmcol removes columns rmvars from table x
% WBC 20200310

rmvars = cellstr(rmvars);
for i = 1:numel(rmvars)
    X.(rmvars{i}) = [];
end

end