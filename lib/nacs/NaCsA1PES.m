function V = NaCsA1PES(r)
    wavenum2hartree = 4.55634011e-6;
    bohr2angstrom = 0.529177210903;
    h = 6.62607003912593e-34;
    hartree = 4.35974464404181e-18;
    
    threshOffs = 4954.24*wavenum2hartree;
    
    A1S_emo = [16501.87,4.65,10509.810,5992.060,...
    4.65431678872284,0.409201927134440,-0.0316930094322587,...
    0.118889237998245,0.272082736776197,0.206826824112031,...
    -0.147778977004416,-0.147778977004416,0.00288158643376113,...
    0.839886790201783,0.358215463426036];
    
    V = wavenum2hartree*emo(r,A1S_emo(1),A1S_emo(2),A1S_emo(3),A1S_emo(4),A1S_emo(5),A1S_emo(6:end)) - threshOffs;

    function U = emo(R,Tdis,rref,Te,De,re,a)
    R = R;
    U = (Tdis - De) + De.*(1 - exp(-alph(R,a,rref).*(R-re))).^2;
    U = U;
    end
    function out = alph(R,a,rref)
        out = 0;
        for i = 1:length(a)
            out = out +  a(i).*((R.^3 - rref.^3)./(R.^3 + rref.^3)).^(i-1);
        end
    end
end

