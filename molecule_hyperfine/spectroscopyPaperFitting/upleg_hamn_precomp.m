%Calculate upleg effective Hamiltonian for different c3Sigma parameters
%with most of the computational heavy lifting already done

function out = upleg_hamn_precomp(pre,c,B)

% build an effective hamiltonian for the interaction between the feshbach
% molecule state, c3sigma, and x1sigma.
if nargin<3
    B = 862e-4;
end

%% load data
%Don't recompute fb unless no file at this field

%% Recompute c3Sigma here
basis = pre.basis;
basis.aUC.ops.Hrot = c.c3Sigma.Be*basis.aUC.ops.J_sq;
basis.aUC.ops.H_OmegaDoubling = c.c3Sigma.wef*pre.basis.aUC.ops.omdub;
basis.aUC.ops.HZ0elecspin = c.c3Sigma.gS*c.uB*pre.basis.aUC.ops.szel;
basis.aUC.ops.Hhf_Na = c.c3Sigma.alpha1*pre.basis.aUC.ops.iNa;
basis.aUC.ops.Hhf_Cs = c.c3Sigma.alpha2*pre.basis.aUC.ops.iCs;
basis.aUC.ops.H = basis.aUC.ops.Hrot ...
    + basis.aUC.ops.HZ0elecspin.*B + basis.aUC.ops.Hhf_Na ...
    + basis.aUC.ops.Hhf_Cs + basis.aUC.ops.H_OmegaDoubling;

H = pagemtimes(basis.change.aUC_FC,'ctranspose',pagemtimes(basis.aUC.ops.H,basis.change.aUC_FC),'none');

pre.c3Sigma.ops.H = H;
[c_psi,evals] = eig(H);
E = real(diag(evals));
out.E = E + c.E_vib;

if 0
    %%
    mF_expect = real(diag(c_psi'*basis.aFC.ops.F_z*c_psi));
    mJ_expect = real(diag(c_psi'*basis.aFC.ops.J_z*c_psi));
    mI_expect = real(diag(c_psi'*basis.aFC.ops.I_z*c_psi));
    figure(12);
    clf;
    hold on; box on;
    for i = 1:numel(E)
        plot(mF_expect(i) + [-0.4 0.4],[1 1]*E(i)*1e-9/c.h,'-k')
    end
    hold off;
    xlabel('M_{tot}');
    ylabel('Energy (GHz)');
end

% rotational TDMs in basis of eigenstates
rot_TDM_f_z_cbasis = 0*pagemtimes(c_psi,'ctranspose',pre.rot_TDM_f_z,'none');
rot_TDM_f_perp_cbasis = pagemtimes(c_psi,'ctranspose',pre.rot_TDM_f_perp,'none');

% up leg hamiltonian
H_up = [];
H_up(:,pre.f_X_rows,:,:) = permute(pre.vibronic_TDM_B_f,[2 1 3 4]).*rot_TDM_f_perp_cbasis(:,pre.f_X_rows,:,:);
H_up(:,pre.f_a_rows,:,:) = permute(pre.vibronic_TDM_b_f,[2 1 3 4]).*rot_TDM_f_perp_cbasis(:,pre.f_a_rows,:,:) ...
    + permute(pre.vibronic_TDM_c_f,[2 1 3 4]).*rot_TDM_f_z_cbasis(:,pre.f_a_rows,:,:);
H_up = permute(sum(H_up,2),[1 3 2 4]);

power = c.power; %W
waist = c.waist; %m

if ~isnumeric(c.pol_up)
    p_up = sphten(polSwitch(c.pol_up));
else
    p_up = sphten(polFrac(c.pol_up));
end

Efield = sqrt(4*c.eta0*power/(pi*waist^2));
out.freq_offs = min(min(E - pre.f.E))/c.h;
out.rabi = Efield*sum(H_up.*pre.psi_init.*permute(p_up,[1 3 4 2]),4);
  
%% save
out.H_up = H_up;
out.B = B;

end

function pOut = polSwitch(pol)
    switch pol
        case 'sigp'
            pOut = [1 -1i 0]/sqrt(2);
        case 'sigm'
            pOut = [1 1i 0]/sqrt(2);
        case 'pi'
            pOut = [0 0 1];
        case 'sigpm'
            pOut = [1 0 0];
    end
end

function pOut = polFrac(sigp_frac)
    sigp_frac = min(sigp_frac,1);
    sigp_frac = max(sigp_frac,0);
    pOut = sqrt(sigp_frac)*[1 -1i 0]/sqrt(2) + sqrt(1 - sigp_frac)*[1 1i 0]/sqrt(2);
end