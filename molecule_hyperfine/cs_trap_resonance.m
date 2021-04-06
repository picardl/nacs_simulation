
% atomic zeeman shifts for Na+Cs

clear;

c = constants();

Cs = atom_basis('Cs');

Cs = couple_angmom(Cs,'s_Cs','i_Cs','f_Cs');

B1 = linspace(0.001e-4,863.7e-4,1e3);

H = Cs.ops.H_Cs_hyperfine+(Cs.ops.H0_Cs_zeeman).*reshape(B1,1,1,[]);

[V,D] = eigenshuffle(H);

f_aobd = [7:7:70];

figure(1); clf;
hold on; box on;
for i = 1:numel(f_aobd)
    plot(B1*1e4,D*1e-6/c.h + f_aobd(i));
end
hold off;
set(gca,'fontsize',14);
xlabel('B (Gauss)');
ylabel('E (MHz)');

% [~,i1] = evec_ind({'f_Cs','m_f_Cs','f_Na','m_f_Na'},[3,3,1,1],NaCs,V(:,:,1));
% [~,i2] = evec_ind({'f_Cs','m_f_Cs','f_Na','m_f_Na'},[3,3,1,1],NaCs,V(:,:,end));

% (D(i2,end) - D(i1,1))/c.h*1e-9

% D(i1,1)/c.h*1e-9

% [~,i1] = evec_ind({'f_Cs','m_f_Cs'},[4,4],Cs,V(:,:,1));
% [~,i2] = evec_ind({'f_Cs','m_f_Cs'},[3,3],Cs,V(:,:,1));
% [~,i3] = evec_ind({'f_Cs','m_f_Cs'},[4,4],Cs,V(:,:,end));
% [~,i4] = evec_ind({'f_Cs','m_f_Cs'},[3,3],Cs,V(:,:,end));
% diff(D([i4 i3],end)/c.h)*1e-9 - diff(D([i2 i1],1)/c.h)*1e-9
% diff(D([i3 i4],1)/c.h)
% (D(i1,1)/c.h - D(i2,1)/c.h) - (D(i3,end)/c.h - D(i4,end)/c.h)

% ((D(evec_ind({'f_Cs','m_f_Cs'},[4,4],Cs,V(:,:,end)),end) - D(evec_ind({'f_Cs','m_f_Cs'},[3,3],Cs,V(:,:,end)),end)) - ...
%     (D(evec_ind({'f_Cs','m_f_Cs'},[4,4],Cs,V(:,:,1)),1) - D(evec_ind({'f_Cs','m_f_Cs'},[3,3],Cs,V(:,:,1)),1)))/c.h

