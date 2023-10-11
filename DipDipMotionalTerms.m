nis = 0:1:4;
njs = 0:1:4;

[nis2D, njs2D] = meshgrid(nis,njs);

xi = 40e-9;
xj = 40e-9;
R0 = 1.6e-6;

[allInt,allResInt] = secondOrder1DTerms(nis2D,njs2D,xi,xj,R0);

normAllInt = allInt./(allInt(1,1))
normAllResInt = allResInt./(allResInt(1,1))

outTab = padarray(normAllResInt,[1 1],0,'pre');
outTab(1,2:end) = nis;
outTab(2:end,1) = nis;

T = array2table(outTab);
table2latex(T)

function [intTotal,intResonant] = secondOrder1DTerms(ni,nj,xi,xj,R0)
    %Zeroth order term
    ch0 = 1/R0^3;

    %First order terms
    ch1i = 3*(ni).*xi/R0^4;
    ch1idag = 3*(ni + 1).*xi/R0^4;
    ch1j = 3*(nj).*xj/R0^4;
    ch1jdag = 3*(nj + 1).*xj/R0^4;

    %Second order terms
    ch2iidag = 6*(ni+1).*xi.^2/R0^5;
    ch2jjdag = 6*(nj+1).*xj.^2/R0^5;
    ch2idagi = 6*(ni).*xi.^2/R0^5;
    ch2jdagj = 6*(nj).*xj.^2/R0^5;
    ch2ij = 2*6*sqrt(ni+1).*sqrt(nj).*xi.*xj/R0^5; %Factor of two for a_ia_j^dag term and h.c.
    ch2ji = 2*6*sqrt(nj+1).*sqrt(ni).*xi.*xj/R0^5;

    intTotal = ch0 - ch1i - ch1idag - ch1j - ch1jdag + ch2iidag + ch2jjdag + ch2idagi + ch2jdagj + ch2ij + ch2ji;
    intResonant = ch0 + ch2iidag + ch2jjdag + ch2idagi + ch2jdagj + ch2ij + ch2ji;
end