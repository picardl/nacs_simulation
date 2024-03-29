function out = wigD_element(r,k,p,q,c,Jvar,Ovar,Mvar)

r_vars = r.Properties.VariableNames;
c_vars = c.Properties.VariableNames;
if ~all(ismember(r_vars,c_vars))
    error('row and column matrix elements must have the same variables');
end

jp = r.(Jvar);
op = r.(Ovar);
mp = r.(Mvar);
j = c.(Jvar);
o = c.(Ovar);
m = c.(Mvar);

r_other = rmcol(rmcol(rmcol(r,Jvar),Mvar),Ovar);
c_other = rmcol(rmcol(rmcol(c,Jvar),Mvar),Ovar);
boo = all(r_other{:,:} == c_other{:,:},2);

out = boo.*(-1).^(mp-op).*sqrt((2*jp+1).*(2*j+1)).*w3j(jp,-op,k,q,j,o).*w3j(jp,-mp,k,p,j,m);

end