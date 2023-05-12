function qnums_cpl = couple_qnums(qnums_unc,j1_name,j2_name,j3_name)

m1_name = ['m_' j1_name];
m2_name = ['m_' j2_name];
m3_name = ['m_' j3_name];

% identify spectator quantum numbers
spec_qnum_names = qnums_unc.Properties.VariableNames(...
    ~ismember(qnums_unc.Properties.VariableNames,...
    {j1_name m1_name j2_name m2_name}));

all_qnum_names = [j1_name j2_name j3_name m3_name spec_qnum_names];

% angular momenta to be coupled
uj1 = unique(qnums_unc.(j1_name));
uj2 = unique(qnums_unc.(j2_name));

% generate values of coupled quantum number
coupled_qnums = zeros(0,size(all_qnum_names,2));
for q = 1:numel(uj1)
    for r = 1:numel(uj2)
        spec_qnums_qr = unique(qnums_unc{qnums_unc.(j1_name)==uj1(q) & qnums_unc.(j2_name)==uj2(r),spec_qnum_names},'rows');
        
        j3new = (abs(uj1(q)-uj2(r)):(uj1(q)+uj2(r)))';
        for s = 1:numel(j3new)
            m3 = (-j3new(s):j3new(s))'; 
            j3 = j3new(s)*ones(size(m3));
            j1 = uj1(q)*ones(size(m3));
            j2 = uj2(r)*ones(size(m3));
            
            coupled_qnums_qr = [j1 j2 j3 m3];
            
            [q_new_ind,q_spec_ind] = ndgrid(1:size(coupled_qnums_qr,1),1:size(spec_qnums_qr,1));
            q_new_ind = q_new_ind(:);
            q_spec_ind = q_spec_ind(:);
            
            coupled_qnums_qr = cat(2,coupled_qnums_qr(q_new_ind,:),spec_qnums_qr(q_spec_ind,:));
            
            coupled_qnums = cat(1,coupled_qnums,coupled_qnums_qr);
        end
    end
end

qnums_cpl = array2table(coupled_qnums,'VariableNames',all_qnum_names);

end