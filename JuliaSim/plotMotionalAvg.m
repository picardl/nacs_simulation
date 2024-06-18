dataDir = 'C:\projects\nacs_simulation\JuliaSim\data\All_n_20231025\';
testFile =  csvread([dataDir,'RamseyData_nxyz1_0_0_0_nxyz2_0_0_0__20231024_171024.csv'],1,0);
% dataDir = 'C:\projects\nacs_simulation\JuliaSim\data\All_n_20231027\';
% testFile =  csvread([dataDir,'RamseyData_nxyz1_0_0_0_nxyz2_0_0_0__20231027_113752.csv'],1,0);
% dataDir = 'C:\projects\nacs_simulation\JuliaSim\data\All_n_20231108\';
% testFile =  csvread([dataDir,'RamseyData_nxyz1_0_0_0_nxyz2_0_0_0__20231108_120759.csv'],1,0);
times = testFile(:,1);

% ground1D = [(1/1.18),(1/1.15),(1/1.13)]; %Based on Na temps in array paper
% ground1D = [1,1,0.9]; %Based on Na temps in array paper
% ground1D = [(1/1.18),(1/1.15),(1/1.13)]*0.5; %Half the expected ground state fraction in each axis
ground1D = [1,1,0.4]; %Based on Na temps in array paper

alpha1D = 1 - ground1D;
ground3D = prod(ground1D);

nx1 = 0:1:2;
ny1 = 0:1:2;
nz1 = 0:1:2;
nx2 = 0:1:2;
ny2 = 0:1:2;
nz2 = 0:1:2;

allP = zeros(length(nx1),length(ny1),length(nz1),length(nx2),length(ny2),length(nz2));
allOsc = zeros(length(nx1),length(ny1),length(nz1),length(nx2),length(ny2),length(nz2),size(testFile,1));
avgOsc = zeros(size(testFile,1),1);

for i1 = 1:length(nx1)
    for i2 =1:length(ny1)
        for i3 = 1:length(nz1)
            for i4 = 1:length(nx2)
                for i5 = 1:length(ny2)
                    for i6 = 1:length(nz2)
                        allP(i1,i2,i3,i4,i5,i6) = probn(nx1(i1),alpha1D(1))*probn(ny1(i2),alpha1D(2))*probn(nz1(i3),alpha1D(3))*...
                            probn(nx2(i4),alpha1D(1))*probn(ny2(i5),alpha1D(2))*probn(nz2(i6),alpha1D(3));
                    end
                end
            end
        end
    end
end

normFac = sum(allP,"all");

for i1 = 1:length(nx1)
    for i2 =1:length(ny1)
        for i3 = 1:length(nz1)
            for i4 = 1:length(nx2)
                for i5 = 1:length(ny2)
                    for i6 = 1:length(nz2)
                        names = dir([dataDir,'RamseyData_nxyz1_',num2str(nx1(i1)),'_',num2str(ny1(i2)),...
                            '_',num2str(nz1(i3)),'_nxyz2_',num2str(nx2(i4)),'_',num2str(ny2(i5)),'_',num2str(nz2(i6)),'__*.csv']);
                        thisOsc =  csvread([dataDir,names.name],1,0);
                        allOsc(i1,i2,i3,i4,i5,i6,:) = thisOsc(:,2);
                        avgOsc = avgOsc + allP(i1,i2,i3,i4,i5,i6)*thisOsc(:,2)/normFac;
                    end
                end
            end
        end
    end
end

figure(999)
% hold on
plot(times,avgOsc)
xlabel({'Interaction time (ms)',dataDir})
ylabel('Joint |00\rangle probability')
title({['R1 p_0 = ',num2str(ground1D(1),2)],['R2 p_0 = ',num2str(ground1D(2),2)],['Ax p_0 = ',num2str(ground1D(3),2)]})

function p = probn(n,alpha)
    p = (1 - alpha)*alpha^n;
end