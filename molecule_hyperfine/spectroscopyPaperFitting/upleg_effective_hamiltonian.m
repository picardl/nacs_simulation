function out = upleg_effective_hamiltonian(const,B,basis,recompute,Jmax,mtot,f_vib)

% build an effective hamiltonian for the interaction between the feshbach
% molecule state, c3sigma, and x1sigma.

if nargin<1
    B = 853e-4;
end
if nargin<2
    basis = 'aFC';
end
if nargin<3
    recompute = 0;
end
if nargin<4
    Jmax = 3;
end
if nargin<5
    mtot = [3 4 5]; 
end
if nargin<6
    f_vib = 325300e9; %Frequency closest to target v level, default v=26
end

%% check for file at this B field
files = dir('../data');
file_ind = contains({files.name},['raman_' strrep(num2str(B*1e4),'.','p') 'G' '_' basis]);
if any(file_ind) && ~(recompute>0)
    disp('found raman hamiltonian file for this B field and basis')
    fnames = {files(file_ind).name};
    times = datenum(regexp(fnames,'\d{6}_\d{6}','match','once'),'YYmmDD_HHMMSS');
    data = load(['../data/' fnames{times==max(times)}]);
    out = data.out;
    return
end

%% load data
f = feshbach(B,basis,1,const);
c = c3Sigma(B,basis,recompute,Jmax,mtot,f_vib,const);

%% ignore quantum numbers that are not common to all 3 states
common_qnums = intersect(c.qnums.Properties.VariableNames,f.qnums.Properties.VariableNames);

% ignore qnums not common to all bases
f.qnums(:,~ismember(f.qnums.Properties.VariableNames,common_qnums)) = [];
c.qnums(:,~ismember(c.qnums.Properties.VariableNames,common_qnums)) = [];

% ignore hund's case (a) quantum numbers in rotational matrix element --
% use only hund's case (c)
qnums_ignore = {'S','Sigma','Lambda'};
f_ignore = ismember(f.qnums.Properties.VariableNames,qnums_ignore);
c_ignore = ismember(c.qnums.Properties.VariableNames,qnums_ignore);


%% transition dipole moments, angular momentum part
p = -1:1; % spherical index

T1q_f_z = sphten([0 0 1]);
T1q_f_perp = sphten([1 1 0]/sqrt(2));

% transition dipole moments
switch basis
    case 'aFC'
        rot_TDM_f_z_components = operator_matrix(@transition_dipole_case_aFC,...
            {c.qnums(:,~c_ignore),f.qnums(:,~f_ignore)},{'eta','J','Omega','I','F','m_F'},p,T1q_f_z);
        rot_TDM_f_perp_components = operator_matrix(@transition_dipole_case_aFC,...
            {c.qnums(:,~c_ignore),f.qnums(:,~f_ignore)},{'eta','J','Omega','I','F','m_F'},p,T1q_f_perp);
    case {'aUC','aIC'}
        rot_TDM_f_z_components = operator_matrix(@transition_dipole_case_a,...
            {c.qnums(:,~c_ignore),f.qnums(:,~f_ignore)},{'eta','J','Omega','m_J'},p,T1q_f_z);
        rot_TDM_f_perp_components = operator_matrix(@transition_dipole_case_a,...
            {c.qnums(:,~c_ignore),f.qnums(:,~f_ignore)},{'eta','J','Omega','m_J'},p,T1q_f_perp);
end

% off-diagonal Hamiltonian blocks
rot_TDM_f_z = (-1).^reshape(p,1,1,3) .* rot_TDM_f_z_components;
rot_TDM_f_perp = (-1).^reshape(p,1,1,3) .* rot_TDM_f_perp_components;

% push polarization info to 4th dimension
rot_TDM_f_z = permute(rot_TDM_f_z,[1,2,4,3,5]);
rot_TDM_f_perp = permute(rot_TDM_f_perp,[1,2,4,3,5]);

%% transition dipole moments, vibronic part
a3S_b3P_elecTDM_data = importdata('lib/rosario_potentials/DM_a3S_b3P_Full');
a3S_c3S_elecTDM_data = importdata('lib/rosario_potentials/DM_a3S_c3S_Full');
X1S_B1P_elecTDM_data = importdata('lib/rosario_potentials/DM_X1S_B1P_Full');

X1S_B1P_elecTDM_interp2 = interp1(X1S_B1P_elecTDM_data(:,1)*const.abohr,X1S_B1P_elecTDM_data(:,2)*const.e*const.abohr,c.r,'spline');
a3S_c3S_elecTDM_interp = interp1(a3S_c3S_elecTDM_data(:,1)*const.abohr,a3S_c3S_elecTDM_data(:,2)*const.e*const.abohr,c.r,'spline');
a3S_b3P_elecTDM_interp = interp1(a3S_b3P_elecTDM_data(:,1)*const.abohr,a3S_b3P_elecTDM_data(:,2)*const.e*const.abohr,c.r,'spline');

%% laser hamiltonian matrix elements

% interpolate radial coord
psi_feshbach_interp = permute(interp1(f.r,permute(f.psi,[2 1 3]),c.r,'spline'),[2 1 3]);

f_X_rows = all(f.qnums{:,{'Lambda','S'}}==[0 0],2);
f_a_rows = all(f.qnums{:,{'Lambda','S'}}==[0 1],2);
c_B_rows = strcmp(cellstr(c.qnums_vib.term),'B1P1');
c_b_rows = strcmp(cellstr(c.qnums_vib.term),'b3P1');
c_c_rows = strcmp(cellstr(c.qnums_vib.term),'c3S1');

vibronic_TDM_B_f = trapz(c.r,psi_feshbach_interp(f_X_rows,:,:).*X1S_B1P_elecTDM_interp2.*c.psi_vib(c_B_rows,:),2);
vibronic_TDM_b_f = trapz(c.r,psi_feshbach_interp(f_a_rows,:,:).*a3S_b3P_elecTDM_interp.*c.psi_vib(c_b_rows,:),2);
vibronic_TDM_c_f = trapz(c.r,psi_feshbach_interp(f_a_rows,:,:).*a3S_c3S_elecTDM_interp.*c.psi_vib(c_c_rows,:),2);

% rotational TDMs in basis of eigenstates
rot_TDM_f_z_cbasis = pagemtimes(c.psi,'ctranspose',rot_TDM_f_z,'none');
rot_TDM_f_perp_cbasis = pagemtimes(c.psi,'ctranspose',rot_TDM_f_perp,'none');

% up leg hamiltonian
H_up(:,f_X_rows,:,:) = permute(vibronic_TDM_B_f,[2 1 3 4]).*rot_TDM_f_perp_cbasis(:,f_X_rows,:,:);
H_up(:,f_a_rows,:,:) = permute(vibronic_TDM_b_f,[2 1 3 4]).*rot_TDM_f_perp_cbasis(:,f_a_rows,:,:) ...
    + permute(vibronic_TDM_c_f,[2 1 3 4]).*rot_TDM_f_z_cbasis(:,f_a_rows,:,:);
H_up = permute(sum(H_up,2),[1 3 2 4]);

%% save
out.f = f;
out.c = c;
out.H_up = H_up;
out.B = B;

if 0
    fn = ['../data/raman_' [strrep(num2str(B*1e4),'.','p') 'G'] '_' basis '_' datestr(now,'YYmmDD_HHMMSS') '.mat'];
    % fn = ['../data/raman_' datestr(now,'YYmmDD_HHMMSS') '.mat'];
    save(fn,'out');
    disp(fn);
end

end
