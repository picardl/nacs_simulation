const = constants();

BFields = [840:5:860,862,864,866,870]*1e-4;
basis = 'aFC';
recompute = 2;
Egap = zeros(length(BFields),1);
rExp = zeros(length(BFields),1);
Emin = zeros(length(BFields),1);

Ethresh =  zeros(length(BFields),1);
lbRat = zeros(length(BFields),1);
E_all = cell(size(BFields));
psi_all = cell(size(BFields));

for i = 1:length(BFields)
    B = BFields(i);
    %% check for file at this B field
    files = dir('../data');
    file_ind = contains({files.name},['fb_' strrep(num2str(B*1e4),'.','p') 'G' '_' basis]);
    if any(file_ind) && ~(recompute>0)
        disp('found fb hamiltonian file for this B field and basis')
        fnames = {files(file_ind).name};
        times = datenum(regexp(fnames,'\d{6}_\d{6}','match','once'),'YYmmDD_HHMMSS');
        data = load(['../data/' fnames{times==max(times)}]);
        out = data.out;
    else
        out = feshbach(B,basis,recompute);
    end
    psi_f = sum(abs(out.psi(:,:,end)),1);
    rExp(i) = trapz(out.r(1:400),out.r(1:400).*psi_f(:,1:400));
    Egap(i) = out.E(end) - out.E(end-1);
    Emin(i) = min(out.E);
    E_all{i} = out.E;
    psi_all{i} = out.psi;
    Ethresh(i) = out.E_lowest_chan_threshold;
end
figure(4)
clf;
plot(BFields*1e4,1./((Emin/const.h - Ethresh/const.h)/1e9),'Marker','x');
xlabel('B Field [G]')
ylabel('1/(E_b - E_{thresh})')
% ylim([-1e2,1e2])
hold on
% plot(BFields*1e4,Emin);
% ylabel('Expectation val of r')
% hold on
% yyaxis right
% scatter(BFields*1e4,Egap);
% ylabel('Gap between atom and fb state')
% xlabel('B Field [G]')

if 0
    %%
    allEs = zeros(8,9);
   for i = 1:length(E_all) 
       allEs(:,i) = E_all{i}(end-7:end) - Ethresh(i)
   end
   figure(4)
   plot(BFields*1e4,allEs)
end