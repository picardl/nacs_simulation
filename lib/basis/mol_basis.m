function q = mol_basis(hunds_case,Lambda_Omega,Nmax_Jmax,S,Omega_fix)

switch hunds_case
    case 'a'
        Lambda = unique([-abs(Lambda_Omega) abs(Lambda_Omega)]);
        Sigma = -max(S):max(S);
        J = 0:Nmax_Jmax;
        m_J = -max(J):max(J);
        [m_J,J,Lambda,S,Sigma] = ndgrid(m_J,J,Lambda,S,Sigma);
        
        q.Lambda = Lambda(:);
        q.S = S(:);
        q.Sigma = Sigma(:);
        q.Omega = q.Lambda+q.Sigma;
        q.J = J(:);
        q.m_J = m_J(:);
        
        q = struct2table(q);
        
        del = abs(q.m_J)>q.J | q.J<abs(q.Omega) | q.S<abs(q.Sigma);
        q(del,:) = [];
        
        if exist('Omega_fix','var')
            q(abs(q.Omega)~=Omega_fix,:) = [];
        end
        
    case 'b'
        Lambda = unique([-abs(Lambda_Omega) abs(Lambda_Omega)]);
        N = abs(Lambda_Omega):Nmax_Jmax;
        m_N = -max(N):max(N);
        
        if exist('S','var')
            m_S = -max(S):max(S);
            [S,m_S,m_N,N,Lambda] = ndgrid(S,m_S,m_N,N,Lambda);
            q.Lambda = Lambda(:);
            q.N = N(:);
            q.m_N = m_N(:);
            q.S = S(:);
            q.m_S = m_S(:);
            
            q = struct2table(q);
            
            del = abs(q.m_N)>q.N | q.N<abs(q.Lambda) | abs(q.m_S)>q.S;
            q(del,:) = [];
            
            q = couple_qnums(q,'N','S','J');
        else
            [m_N,N,Lambda] = ndgrid(m_N,N,Lambda);
            q.Lambda = Lambda(:);
            q.N = N(:);
            q.m_N = m_N(:);
            
            q = struct2table(q);
            
            del = abs(q.m_N)>q.N | q.N<abs(q.Lambda);
            q(del,:) = [];
        end
        
    case 'c'
        Omega = unique([-abs(Lambda_Omega) abs(Lambda_Omega)]);
        J = abs(Lambda_Omega):Nmax_Jmax;
        m_J = -max(J):max(J);
        [m_J,J,Omega] = ndgrid(m_J,J,Omega);
        q.Omega = Omega(:);
        q.J = J(:);
        q.m_J = m_J(:);
        
        
        q = struct2table(q);
        
        del = abs(q.m_J)>q.J | q.J<abs(q.Omega);
        q(del,:) = [];
        
    otherwise
        error('that hunds case either doesnt exist or hasnt been programmed yet')
end

end