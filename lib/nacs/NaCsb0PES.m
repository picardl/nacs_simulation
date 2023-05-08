function V = NaCsb0PES(r)
    wavenum2hartree = 4.55634011e-6;
    bohr2angstrom = 0.529177210903;
    h = 6.62607003912593e-34;
    hartree = 4.35974464404181e-18;
    
    threshOffs = 184.68*wavenum2hartree + 4954.24*wavenum2hartree;
    
    b3p_emo = [16317.19,3.78,10236.048,6081.142,...
    3.77999922692013,0.665828316447219,0.159655197110690,...
    0.144813242483051,0.0964993864916689,-0.0452724688520158,...
    0.320383473138789,0.0148834980740537,-0.0180900998067171,...
    -0.388373323943793];
    
    V = wavenum2hartree*emo(r,b3p_emo(1),b3p_emo(2),b3p_emo(3),b3p_emo(4),b3p_emo(5),b3p_emo(6:end)) - threshOffs;

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

