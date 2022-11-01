from lib.am_utils import *
import lib.consts as c

# nacs
p0 = {'I1': 3/2,
      'I2': 7/2,
      'N': arange(0,1),
      'g1': 1.478,
      'g2': 0.748,
      'Bv': 1.7396*1e6, # MHz
      'eQq1': -0.097*1e3, # MHz
      'eQq2': 0.15*1e3, # MHz
      'sig1': 639.2 / 1e6, # ppm
      'sig2': 6278.7 / 1e6, # ppm
      'c1': 14.2*1e-3, # Hz
      'c2': 854.5*1e-3, # Hz
      'c3': 105.6*1e-3, # Hz
      'c4': 3941.8*1e-3, # Hz
      'gr': 0, # 0.0144 ##
      'mu': 4.6*3.33564*1e-28/c.h, # D
      'muN': 5.050783699e-27 * 1e-4 * 1e-3/ c.h, # J/T
      'muE': 4.6*(3.33564e-30), # Debye
      'alpha_p': 1872.1153,
      'alpha_s': 467.0375,
      'E_td_00': np.linspace(0,0,1), # 1000
      'trap_pol': np.array([[0],[1],[0]]),
      'E_dc': np.array([0,0,0]), # V/cm
      'B': np.linspace(866,866,1)}

p0s = {'NaRb': {'I1': 3/2,
          'I2': 3/2,
          'N': arange(0,1),
          'g1': 1.484,
          'g2': 1.832,
          'Bv': 2.0896628*1e6, # MHz
          'eQq1': -0.139*1e3, # MHz
          'eQq2': -3.048*1e3, # .3*1e3, #MHz
          'sig1': 0, # ppm, bundled into g1
          'sig2': 0, # ppm, bundled into g2
          'c1': 60.7*1e-3, # Hz
          'c2': 983.8*1e-3, # Hz
          'c3': 259.3*1e-3, # Hz
          'c4': 6.56, # kHz
          'gr': 0.001, # 0.0144 ##
          'mu': 4.6*3.33564*1e-28/c.h, # D
          'muN': 5.050783699e-27 * 1e-4 * 1e-3/ c.h, # J/T
          'muE': 3.3*(3.33564e-30), # Debye
          'alpha_p': 699.0+469.0,
          'alpha_s': (2*699.0-469.0)/2,
          'E_td_00': np.linspace(0,0,1), # 1000
          'trap_pol': np.array([[0],[0],[1]]),
          'E_dc': np.array([0,0,0]), # V/cm
          'B': np.linspace(0,0,1)}, 
      'NaCs': {'I1': 3/2,
          'I2': 7/2,
          'N': arange(0,1),
          'g1': 1.478,
          'g2': 0.748,
          'Bv': 1.7396*1e6, # MHz
          'eQq1': -0.097*1e3, # MHz
          'eQq2': 0.15*1e3, # MHz
          'sig1': 639.2 / 1e6, # ppm
          'sig2': 6278.7 / 1e6, # ppm
          'c1': 14.2*1e-3, # Hz
          'c2': 854.5*1e-3, # Hz
          'c3': 105.6*1e-3, # Hz
          'c4': 3941.8*1e-3, # Hz
          'gr': 0, # 0.0144 ##
          'mu': 4.6*3.33564*1e-28/c.h, # D
          'muN': 5.050783699e-27 * 1e-4 * 1e-3/ c.h, # J/T
          'muE': 4.6*(3.33564e-30), # Debye
          'alpha_p': 1872.1153,
          'alpha_s': 467.0375,
          'E_td_00': np.linspace(0,0,1), # 1000
          'trap_pol': np.array([[0],[1],[0]]),
          'E_dc': np.array([0,0,0]), # V/cm
          'B': np.linspace(0,0,1)},
      'KRb': {'I1': 4,
          'I2': 3/2,
          'N': arange(0,1),
          'g1': -0.324,
          'g2': 1.834,
          'Bv': 1.114*1e6, # MHz
          'eQq1': 0.306*1e3, # MHz
          'eQq2': -1.520*1e3, # MHz
          'sig1': 1321 / 1e6, # ppm
          'sig2': 3469 / 1e6, # ppm
          'c1': -24.1*1e-3, # Hz
          'c2': 420.1*1e-3, # Hz
          'c3': -48.2*1e-3, # Hz
          'c4': 2030.4*1e-3, # Hz
          'gr': 0.0144, # 0.0144 ##
          'mu': 4.6*3.33564*1e-28/c.h, # D
          'muN': 5.050783699e-27 * 1e-4 * 1e-3/ c.h, # J/T
          'muE': 0.76*(3.33564e-30), # Debye
          'alpha_p': 10, # in intensity units, but factors out so only ratio matters
          'alpha_s': 3.3,
          'E_td_00': np.linspace(0,0,1), # 1000
          'trap_pol': np.array([[0],[1],[0]]),
          'E_dc': np.array([0,0,0]), # V/cm
          'B': np.linspace(0,0,1)}, 
      'RbCs': {'I1': 3/2,
          'I2': 7/2,
          'N': arange(0,1),
          'g1': 1.834,
          'g2': 0.738,
          'Bv': 0.490173*1e6, # MHz
          'eQq1': -0.80929*1e3, # MHz
          'eQq2': 0.05998*1e3, # MHz
          'sig1': 3531 / 1e6, # ppm
          'sig2': 6367 / 1e6, # ppm
          'c1': 29.4*1e-3, # Hz
          'c2': 194.1*1e-3, # Hz
          'c3': 192.4*1e-3, # Hz
          'c4': 19.019, # Hz
          'gr': 0.0062, # 0.0144 ##
          'muE': 1.25*3.33564*1e-28/c.h / 1e3, # D
          'muN': 5.050783699e-27 * 1e-4 * 1e-3/ c.h, # J/T
          'alpha_p': 2020+1997, # 1997*3/2, # 1872.1153,
          'alpha_s': (2*2020-1997)/2, # 0, # 467.0375,
          'E_td_00': np.linspace(0,0,1), # 1000
          'trap_pol': np.array([[0],[1],[0]]),
          'E_dc': np.array([0,0,0]), # V/cm
          'B': np.linspace(0,0,1)}}

def dm_ops_sph(b,j='n'):
    _,_,delta_spec = get_spectators(b,[j,'m_'+j],mat=True)
    
    f1 = b_op(lambda b1,b2,j: (-1)**(b1['m_'+j])*np.sqrt((2*b1[j]+1)*(2*b2[j]+1)), b, 0,j) 
    f2 = b_op(lambda b1,b2,j: py3nj.wigner3j(b1[j], int(2*1), b2[j],
                                             0,     0,        0), b,1,j)
    
    cg_z = b_op(lambda b1,b2,j,q: py3nj.wigner3j(b1[j],    int(2*1),   b2[j],
                                              -b1['m_'+j], int(2*q), b2['m_'+j]), b,1,j,0)
    cg_r1 = b_op(lambda b1,b2,j,q: py3nj.wigner3j(b1[j],   int(2*1),   b2[j],
                                              -b1['m_'+j], int(2*q), b2['m_'+j]), b,1,j,1)
    cg_r2 = b_op(lambda b1,b2,j,q: py3nj.wigner3j(b1[j],   int(2*1),   b2[j],
                                              -b1['m_'+j], int(2*q), b2['m_'+j]), b,1,j,-1)
        
    M_z = delta_spec*f1*f2*cg_z
    M_r1 = delta_spec*f1*f2*cg_r1 # sig+
    M_r2 = delta_spec*f1*f2*cg_r2 # sig-

    return (M_r2,M_z,M_r1)

def H_mol(p):
    
    # create basis
    b_uc = am_basis([['i1',p['I1'],1],['i2',p['I2'],1],['n',p['N'],1]])
    b_sc, M_uc2sc = couple_basis(b_uc,['i1','i2','i'])
    b_fic, M_uc2fic = ([], [])
    for i,f in zip(['i1','i2'],['f1','f2']):
        b_t, M_t = couple_basis(b_uc,[i,'n',f])
        b_fic.append(b_t)
        M_uc2fic.append(M_t)
    b_fc, M_sc2fc = couple_basis(b_sc,['i','n','f'])
    M_uc2fc = M_sc2fc@M_uc2sc    
    

    H_hf_eQi_fic, H_hf_eQi_fc, H_hf_eQi_uc = ([],[],[])
    for ind,j in enumerate(zip(['i1','i2'],['f1','f2'])):
        i,n,f = (j[0],'n',j[1])
        v = b_fic[ind].v.copy()
        C = np.diag(op(b_fic[ind],i+'*'+n+','+f))
        H_hf_eQi_fic.append(np.diag(((3*C*C + 3/2*C - v[i]*(v[i]+1)*v[n]*(v[n]+1)) /
                      (2*v[i]*(2*v[i]-1)*(2*v[n]-1)*(2*v[n]+3)))))
        #H_hf_eQi_fic.append(np.diag(((3/2*C*(2*C+1) - v[i]*(v[i]+1)*v[n]*(v[n]+1)) /
        #              (2*v[i]*(2*v[i]-1)*(2*v[n]-1)*(2*v[n]+3)))))
        H_hf_eQi_uc.append(np.transpose(M_uc2fic[ind])@H_hf_eQi_fic[ind]@M_uc2fic[ind])
        H_hf_eQi_fc.append(M_uc2fc@np.transpose(M_uc2fic[ind])@H_hf_eQi_fic[ind]@
                             M_uc2fic[ind]@np.transpose(M_uc2fc))
    H_hf_eQ_uc = - p['eQq1']*H_hf_eQi_uc[0] - p['eQq2']*H_hf_eQi_uc[1]

    H_rot_uc = p['Bv'] * op(b_uc,'n^2')
    H_z_rot_uc = -p['gr']*p['muN'] * op(b_uc,'m_n')
    H_z_nuc_uc = -p['g1']*p['muN']*(1-p['sig1']) * op(b_uc,'m_i1') \
                    -p['g2']*p['muN']*(1-p['sig2']) * op(b_uc,'m_i2')

    H_hf_sr_uc = (p['c1'] * np.transpose(M_uc2fic[0])@op(b_fic[0],'i1*n,f1')@M_uc2fic[0] + \
                  p['c2'] * np.transpose(M_uc2fic[1])@op(b_fic[1],'i2*n,f2')@M_uc2fic[1])
    H_hf_ss_uc = p['c4'] * np.transpose(M_uc2fc)@op(b_fc,'i1*i2,i')@M_uc2fc
    
#     def H_ls(b,a_p,a_s,E00,pol):
        
#         def rot_state_op(b,M,a00):
#             Mrot = b_op(lambda b1,b2: ((b1['n']==0)&(b2['n']==0))*a00 + 
#             (((b1['n']==1)&(b2['n']==1)&(b1['m_n']==-1)&(b2['m_n']==-1))*M[0,0] +
#              ((b1['n']==1)&(b2['n']==1)&(b1['m_n']==-1)&(b2['m_n']==0))*M[0,1] +
#              ((b1['n']==1)&(b2['n']==1)&(b1['m_n']==-1)&(b2['m_n']==1))*M[0,2] +
#              ((b1['n']==1)&(b2['n']==1)&(b1['m_n']==0)&(b2['m_n']==-1))*M[1,0] +
#              ((b1['n']==1)&(b2['n']==1)&(b1['m_n']==0)&(b2['m_n']==0))*M[1,1] +
#              ((b1['n']==1)&(b2['n']==1)&(b1['m_n']==0)&(b2['m_n']==1))*M[1,2] +
#              ((b1['n']==1)&(b2['n']==1)&(b1['m_n']==1)&(b2['m_n']==-1))*M[2,0] +
#              ((b1['n']==1)&(b2['n']==1)&(b1['m_n']==1)&(b2['m_n']==0))*M[2,1] +
#              ((b1['n']==1)&(b2['n']==1)&(b1['m_n']==1)&(b2['m_n']==1))*M[2,2]),b,0)
#             return Mrot
        
#         _,_,delta_spec = get_spectators(b,['n','m_n'],mat=True)
#         a_00 = 1/3 * (a_p + 2*a_s)
        
#         a_23 = (2*a_p+3*a_s)/5
#         a_32 = (3*a_p+2*a_s)/5        
#         a_14 = (a_p+4*a_s)/5
#         a_11 = (a_p-a_s)/5
        
#         MM = np.zeros((3,3,3,3),dtype = 'complex_')
#         MM[0,0] = np.array([[a_23,0,0],[0,a_23,0],[0,0,a_14]])
#         MM[0,1] = np.array([[0,0,1j*a_11/np.sqrt(2)],[0,0,a_11/np.sqrt(2)],[1j*a_11/np.sqrt(2),a_11/np.sqrt(2),0]])
#         MM[0,2] = np.array([[a_11,-1j*a_11,0],[-1j*a_11,-a_11,0],[0,0,0]])
#         MM[1,0] = np.array([[0,0,-1j*a_11/np.sqrt(2)],[0,0,a_11/np.sqrt(2)],[-1j*a_11/np.sqrt(2),a_11/np.sqrt(2),0]])
#         MM[1,1] = np.array([[a_14,0,0],[0,a_14,0],[0,0,a_32]])
#         MM[1,2] = np.array([[0,0,-1j*a_11/np.sqrt(2)],[0,0,-a_11/np.sqrt(2)],[-1j*a_11/np.sqrt(2),-a_11/np.sqrt(2),0]])
#         MM[2,0] = np.array([[a_11,1j*a_11,0],[1j*a_11,-a_11,0],[0,0,0]])
#         MM[2,1] = np.array([[0,0,1j*a_11/np.sqrt(2)],[0,0,-a_11/np.sqrt(2)],[1j*a_11/np.sqrt(2),-a_11/np.sqrt(2),0]])
#         MM[2,2] = np.array([[a_23,0,0],[0,a_23,0],[0,0,a_14]])   
#         MM = np.einsum('ijkl->klij',MM)
# #         for i in range(3):
# #             for j in range(3):
# #                 np.fill_diagonal(MM[i,j,:,:],np.diagonal(MM[i,j])-a_00)   
# #                 np.fill_diagonal(MM[i,j,:,:],np.diagonal(MM[i,j])-a_00)        
#         if pol.shape == (3,1):
#             P = np.einsum('ijkl,jm->ikl',MM,pol)
#             Q = np.einsum('ij,jlm->ilm',np.conjugate(pol.T),P)
#         else:
#             P = np.einsum('ijkl,jm->ikl',MM,pol.T)
#             Q = np.einsum('ij,jlm->ilm',np.conjugate(pol),P)
#         Q = Q / np.linalg.norm(pol)
#         mat = rot_state_op(b,Q[0,:],a_00)
        
#         return delta_spec*mat*E00/a_00     

#     H_ls_uc = -H_ls(b_uc,p['alpha_p'],p['alpha_s'],1,p['trap_pol'])

    def H_ls(b,a_p,a_s,E00,pol):    
         
        def rot_state_op(b,M,a00):
            N_mn = M.shape[0]
            if N_mn == 1:
                N_rot = 1
            elif N_mn == 4:
                N_rot = 2
            elif N_mn == 9:
                N_rot = 3
            else:
                print('something wrong')

            rot_state_list = []
            k = 0
            for i in range(N_rot):
                for j in range(-i,i+1):
                    rot_state_list.append([k,i,j])
                    k += 1
            
            def rot_mat_el(b1,b2):
                mat_el = 0
                for r1 in rot_state_list:
                    for r2 in rot_state_list:
                        mat_el += ((b1['n']==r1[1])&(b2['n']==r2[1])&(b1['m_n']==r1[2])&(b2['m_n']==r2[2]))*M[r1[0],r2[0]]
                return mat_el            
            
            Mrot = b_op(rot_mat_el,b,0)
            return Mrot
        
        _,_,delta_spec = get_spectators(b,['n','m_n'],mat=True)
        a_00 = 1/3 * (a_p + 2*a_s)        
        
        if np.max(p['N']) < 2:
            a_23 = (2*a_p+3*a_s)/5
            a_32 = (3*a_p+2*a_s)/5        
            a_14 = (a_p+4*a_s)/5
            a_11 = (a_p-a_s)/5
            
            MM = np.zeros((4,4,3,3),dtype = 'complex_')
            MM[0,0] = np.array([[a_00,0,0],[0,a_00,0],[0,0,a_00]])
            MM[1,1] = np.array([[a_23,0,0],[0,a_23,0],[0,0,a_14]])
            MM[1,2] = np.array([[0,0,1j*a_11/np.sqrt(2)],[0,0,a_11/np.sqrt(2)],[1j*a_11/np.sqrt(2),a_11/np.sqrt(2),0]])
            MM[1,3] = np.array([[a_11,-1j*a_11,0],[-1j*a_11,-a_11,0],[0,0,0]])
            MM[2,1] = np.array([[0,0,-1j*a_11/np.sqrt(2)],[0,0,a_11/np.sqrt(2)],[-1j*a_11/np.sqrt(2),a_11/np.sqrt(2),0]])
            MM[2,2] = np.array([[a_14,0,0],[0,a_14,0],[0,0,a_32]])
            MM[2,3] = np.array([[0,0,-1j*a_11/np.sqrt(2)],[0,0,-a_11/np.sqrt(2)],[-1j*a_11/np.sqrt(2),-a_11/np.sqrt(2),0]])
            MM[3,1] = np.array([[a_11,1j*a_11,0],[1j*a_11,-a_11,0],[0,0,0]])
            MM[3,2] = np.array([[0,0,1j*a_11/np.sqrt(2)],[0,0,-a_11/np.sqrt(2)],[1j*a_11/np.sqrt(2),-a_11/np.sqrt(2),0]])
            MM[3,3] = np.array([[a_23,0,0],[0,a_23,0],[0,0,a_14]])   
            
        elif np.max(p['N']) == 2:
            
            # produced from polarizabilityVsAngle_n2_topython.nb
            MM = np.array( [np.array( [np.array( [np.array( [1/3 * ( a_s + \
2 * a_s ),0,0,] ),np.array( [0,1/3 * ( a_s + 2 * a_s ),0,] \
),np.array( [0,0,1/3 * ( a_s + 2 * a_s ),] ),] ),np.array( \
[np.array( [0,0,0,] ),np.array( [0,0,0,] ),np.array( \
[0,0,0,] ),] ),np.array( [np.array( [0,0,0,] ),np.array( \
[0,0,0,] ),np.array( [0,0,0,] ),] ),np.array( [np.array( \
[0,0,0,] ),np.array( [0,0,0,] ),np.array( [0,0,0,] ),] \
),np.array( [np.array( [( 30 )**( -1/2 ) * ( -1 * a_s + a_s \
),complex( 0,-1 ) * ( 30 )**( -1/2 ) * ( a_s + -1 * a_s ),0,] \
),np.array( [complex( 0,-1 ) * ( 30 )**( -1/2 ) * ( a_s + -1 * a_s \
),( 30 )**( -1/2 ) * ( a_s + -1 * a_s ),0,] ),np.array( [0,0,0,] ),] \
),np.array( [np.array( [0,0,complex( 0,-1 ) * ( 30 )**( -1/2 ) \
* ( a_s + -1 * a_s ),] ),np.array( [0,0,( 30 )**( -1/2 ) * ( a_s + -1 \
* a_s ),] ),np.array( [complex( 0,-1 ) * ( 30 )**( -1/2 ) * ( a_s + \
-1 * a_s ),( 30 )**( -1/2 ) * ( a_s + -1 * a_s ),0,] ),] ),np.array( \
[np.array( [1/3 * ( 5 )**( -1/2 ) * ( -1 * a_s + a_s ),0,0,] \
),np.array( [0,1/3 * ( 5 )**( -1/2 ) * ( -1 * a_s + a_s ),0,] \
),np.array( [0,0,2/3 * ( 5 )**( -1/2 ) * ( a_s + -1 * a_s ),] ),] \
),np.array( [np.array( [0,0,complex( 0,-1 ) * ( 30 )**( -1/2 ) \
* ( a_s + -1 * a_s ),] ),np.array( [0,0,( 30 )**( -1/2 ) * ( -1 * a_s \
+ a_s ),] ),np.array( [complex( 0,-1 ) * ( 30 )**( -1/2 ) * ( a_s + \
-1 * a_s ),( 30 )**( -1/2 ) * ( -1 * a_s + a_s ),0,] ),] ),np.array( \
[np.array( [( 30 )**( -1/2 ) * ( -1 * a_s + a_s ),complex( 0,1 ) * ( \
30 )**( -1/2 ) * ( a_s + -1 * a_s ),0,] ),np.array( [complex( 0,1 ) \
* ( 30 )**( -1/2 ) * ( a_s + -1 * a_s ),( 30 )**( -1/2 ) * ( a_s + -1 * \
a_s ),0,] ),np.array( [0,0,0,] ),] ),] ),np.array( [np.array( \
[np.array( [0,0,0,] ),np.array( [0,0,0,] ),np.array( \
[0,0,0,] ),] ),np.array( [np.array( [1/5 * ( 2 * a_s + 3 * a_s \
),0,0,] ),np.array( [0,1/5 * ( 2 * a_s + 3 * a_s ),0,] \
),np.array( [0,0,1/5 * ( a_s + 4 * a_s ),] ),] ),np.array( \
[np.array( [0,0,complex( 0,1/5 ) * ( 2 )**( -1/2 ) * ( a_s + -1 * \
a_s ),] ),np.array( [0,0,1/5 * ( 2 )**( -1/2 ) * ( a_s + -1 * a_s ),] \
),np.array( [complex( 0,1/5 ) * ( 2 )**( -1/2 ) * ( a_s + -1 * a_s \
),1/5 * ( 2 )**( -1/2 ) * ( a_s + -1 * a_s ),0,] ),] ),np.array( \
[np.array( [1/5 * ( a_s + -1 * a_s ),complex( 0,-1/5 ) * ( a_s + -1 * \
a_s ),0,] ),np.array( [complex( 0,-1/5 ) * ( a_s + -1 * a_s ),1/5 * ( \
-1 * a_s + a_s ),0,] ),np.array( [0,0,0,] ),] ),np.array( \
[np.array( [0,0,0,] ),np.array( [0,0,0,] ),np.array( \
[0,0,0,] ),] ),np.array( [np.array( [0,0,0,] ),np.array( \
[0,0,0,] ),np.array( [0,0,0,] ),] ),np.array( [np.array( \
[0,0,0,] ),np.array( [0,0,0,] ),np.array( [0,0,0,] ),] \
),np.array( [np.array( [0,0,0,] ),np.array( [0,0,0,] \
),np.array( [0,0,0,] ),] ),np.array( [np.array( [0,0,0,] \
),np.array( [0,0,0,] ),np.array( [0,0,0,] ),] ),] \
),np.array( [np.array( [np.array( [0,0,0,] ),np.array( \
[0,0,0,] ),np.array( [0,0,0,] ),] ),np.array( [np.array( \
[0,0,complex( 0,-1/5 ) * ( 2 )**( -1/2 ) * ( a_s + -1 * a_s ),] \
),np.array( [0,0,1/5 * ( 2 )**( -1/2 ) * ( a_s + -1 * a_s ),] \
),np.array( [complex( 0,-1/5 ) * ( 2 )**( -1/2 ) * ( a_s + -1 * a_s \
),1/5 * ( 2 )**( -1/2 ) * ( a_s + -1 * a_s ),0,] ),] ),np.array( \
[np.array( [1/5 * ( a_s + 4 * a_s ),0,0,] ),np.array( [0,1/5 * ( \
a_s + 4 * a_s ),0,] ),np.array( [0,0,1/5 * ( 3 * a_s + 2 * a_s ),] ),] \
),np.array( [np.array( [0,0,complex( 0,-1/5 ) * ( 2 )**( -1/2 ) \
* ( a_s + -1 * a_s ),] ),np.array( [0,0,1/5 * ( 2 )**( -1/2 ) * ( -1 \
* a_s + a_s ),] ),np.array( [complex( 0,-1/5 ) * ( 2 )**( -1/2 ) * ( \
a_s + -1 * a_s ),1/5 * ( 2 )**( -1/2 ) * ( -1 * a_s + a_s ),0,] ),] \
),np.array( [np.array( [0,0,0,] ),np.array( [0,0,0,] \
),np.array( [0,0,0,] ),] ),np.array( [np.array( [0,0,0,] \
),np.array( [0,0,0,] ),np.array( [0,0,0,] ),] ),np.array( \
[np.array( [0,0,0,] ),np.array( [0,0,0,] ),np.array( \
[0,0,0,] ),] ),np.array( [np.array( [0,0,0,] ),np.array( \
[0,0,0,] ),np.array( [0,0,0,] ),] ),np.array( [np.array( \
[0,0,0,] ),np.array( [0,0,0,] ),np.array( [0,0,0,] ),] ),] \
),np.array( [np.array( [np.array( [0,0,0,] ),np.array( \
[0,0,0,] ),np.array( [0,0,0,] ),] ),np.array( [np.array( \
[1/5 * ( a_s + -1 * a_s ),complex( 0,1/5 ) * ( a_s + -1 * a_s ),0,] \
),np.array( [complex( 0,1/5 ) * ( a_s + -1 * a_s ),1/5 * ( -1 * a_s + \
a_s ),0,] ),np.array( [0,0,0,] ),] ),np.array( [np.array( \
[0,0,complex( 0,1/5 ) * ( 2 )**( -1/2 ) * ( a_s + -1 * a_s ),] \
),np.array( [0,0,1/5 * ( 2 )**( -1/2 ) * ( -1 * a_s + a_s ),] \
),np.array( [complex( 0,1/5 ) * ( 2 )**( -1/2 ) * ( a_s + -1 * a_s \
),1/5 * ( 2 )**( -1/2 ) * ( -1 * a_s + a_s ),0,] ),] ),np.array( \
[np.array( [1/5 * ( 2 * a_s + 3 * a_s ),0,0,] ),np.array( [0,1/5 \
* ( 2 * a_s + 3 * a_s ),0,] ),np.array( [0,0,1/5 * ( a_s + 4 * a_s ),] \
),] ),np.array( [np.array( [0,0,0,] ),np.array( [0,0,0,] \
),np.array( [0,0,0,] ),] ),np.array( [np.array( [0,0,0,] \
),np.array( [0,0,0,] ),np.array( [0,0,0,] ),] ),np.array( \
[np.array( [0,0,0,] ),np.array( [0,0,0,] ),np.array( \
[0,0,0,] ),] ),np.array( [np.array( [0,0,0,] ),np.array( \
[0,0,0,] ),np.array( [0,0,0,] ),] ),np.array( [np.array( \
[0,0,0,] ),np.array( [0,0,0,] ),np.array( [0,0,0,] ),] ),] \
),np.array( [np.array( [np.array( [( 30 )**( -1/2 ) * ( -1 * \
a_s + a_s ),complex( 0,1 ) * ( 30 )**( -1/2 ) * ( a_s + -1 * a_s ),0,] \
),np.array( [complex( 0,1 ) * ( 30 )**( -1/2 ) * ( a_s + -1 * a_s \
),( 30 )**( -1/2 ) * ( a_s + -1 * a_s ),0,] ),np.array( [0,0,0,] ),] \
),np.array( [np.array( [0,0,0,] ),np.array( [0,0,0,] \
),np.array( [0,0,0,] ),] ),np.array( [np.array( [0,0,0,] \
),np.array( [0,0,0,] ),np.array( [0,0,0,] ),] ),np.array( \
[np.array( [0,0,0,] ),np.array( [0,0,0,] ),np.array( \
[0,0,0,] ),] ),np.array( [np.array( [1/7 * ( 3 * a_s + 4 * a_s \
),0,0,] ),np.array( [0,1/7 * ( 3 * a_s + 4 * a_s ),0,] \
),np.array( [0,0,1/7 * ( a_s + 6 * a_s ),] ),] ),np.array( \
[np.array( [0,0,complex( 0,1/7 ) * ( a_s + -1 * a_s ),] \
),np.array( [0,0,1/7 * ( a_s + -1 * a_s ),] ),np.array( [complex( \
0,1/7 ) * ( a_s + -1 * a_s ),1/7 * ( a_s + -1 * a_s ),0,] ),] \
),np.array( [np.array( [1/7 * ( 2/3 )**( 1/2 ) * ( a_s + -1 * a_s \
),complex( 0,-1/7 ) * ( 2/3 )**( 1/2 ) * ( a_s + -1 * a_s ),0,] \
),np.array( [complex( 0,-1/7 ) * ( 2/3 )**( 1/2 ) * ( a_s + -1 * a_s \
),1/7 * ( 2/3 )**( 1/2 ) * ( -1 * a_s + a_s ),0,] ),np.array( \
[0,0,0,] ),] ),np.array( [np.array( [0,0,0,] ),np.array( \
[0,0,0,] ),np.array( [0,0,0,] ),] ),np.array( [np.array( \
[0,0,0,] ),np.array( [0,0,0,] ),np.array( [0,0,0,] ),] ),] \
),np.array( [np.array( [np.array( [0,0,complex( 0,1 ) * ( 30 \
)**( -1/2 ) * ( a_s + -1 * a_s ),] ),np.array( [0,0,( 30 )**( -1/2 ) \
* ( a_s + -1 * a_s ),] ),np.array( [complex( 0,1 ) * ( 30 )**( -1/2 \
) * ( a_s + -1 * a_s ),( 30 )**( -1/2 ) * ( a_s + -1 * a_s ),0,] ),] \
),np.array( [np.array( [0,0,0,] ),np.array( [0,0,0,] \
),np.array( [0,0,0,] ),] ),np.array( [np.array( [0,0,0,] \
),np.array( [0,0,0,] ),np.array( [0,0,0,] ),] ),np.array( \
[np.array( [0,0,0,] ),np.array( [0,0,0,] ),np.array( \
[0,0,0,] ),] ),np.array( [np.array( [0,0,complex( 0,-1/7 ) * ( \
a_s + -1 * a_s ),] ),np.array( [0,0,1/7 * ( a_s + -1 * a_s ),] \
),np.array( [complex( 0,-1/7 ) * ( a_s + -1 * a_s ),1/7 * ( a_s + -1 \
* a_s ),0,] ),] ),np.array( [np.array( [1/7 * ( 2 * a_s + 5 * a_s \
),0,0,] ),np.array( [0,1/7 * ( 2 * a_s + 5 * a_s ),0,] \
),np.array( [0,0,1/7 * ( 3 * a_s + 4 * a_s ),] ),] ),np.array( \
[np.array( [0,0,complex( 0,1/7 ) * ( 6 )**( -1/2 ) * ( a_s + -1 * \
a_s ),] ),np.array( [0,0,1/7 * ( 6 )**( -1/2 ) * ( a_s + -1 * a_s ),] \
),np.array( [complex( 0,1/7 ) * ( 6 )**( -1/2 ) * ( a_s + -1 * a_s \
),1/7 * ( 6 )**( -1/2 ) * ( a_s + -1 * a_s ),0,] ),] ),np.array( \
[np.array( [1/7 * ( a_s + -1 * a_s ),complex( 0,-1/7 ) * ( a_s + -1 * \
a_s ),0,] ),np.array( [complex( 0,-1/7 ) * ( a_s + -1 * a_s ),1/7 * ( \
-1 * a_s + a_s ),0,] ),np.array( [0,0,0,] ),] ),np.array( \
[np.array( [0,0,0,] ),np.array( [0,0,0,] ),np.array( \
[0,0,0,] ),] ),] ),np.array( [np.array( [np.array( [1/3 * ( \
5 )**( -1/2 ) * ( -1 * a_s + a_s ),0,0,] ),np.array( [0,1/3 * ( 5 \
)**( -1/2 ) * ( -1 * a_s + a_s ),0,] ),np.array( [0,0,2/3 * ( 5 )**( \
-1/2 ) * ( a_s + -1 * a_s ),] ),] ),np.array( [np.array( [0,0,0,] \
),np.array( [0,0,0,] ),np.array( [0,0,0,] ),] ),np.array( \
[np.array( [0,0,0,] ),np.array( [0,0,0,] ),np.array( \
[0,0,0,] ),] ),np.array( [np.array( [0,0,0,] ),np.array( \
[0,0,0,] ),np.array( [0,0,0,] ),] ),np.array( [np.array( \
[1/7 * ( 2/3 )**( 1/2 ) * ( a_s + -1 * a_s ),complex( 0,1/7 ) * ( 2/3 \
)**( 1/2 ) * ( a_s + -1 * a_s ),0,] ),np.array( [complex( 0,1/7 ) * \
( 2/3 )**( 1/2 ) * ( a_s + -1 * a_s ),1/7 * ( 2/3 )**( 1/2 ) * ( -1 * \
a_s + a_s ),0,] ),np.array( [0,0,0,] ),] ),np.array( \
[np.array( [0,0,complex( 0,-1/7 ) * ( 6 )**( -1/2 ) * ( a_s + -1 * \
a_s ),] ),np.array( [0,0,1/7 * ( 6 )**( -1/2 ) * ( a_s + -1 * a_s ),] \
),np.array( [complex( 0,-1/7 ) * ( 6 )**( -1/2 ) * ( a_s + -1 * a_s \
),1/7 * ( 6 )**( -1/2 ) * ( a_s + -1 * a_s ),0,] ),] ),np.array( \
[np.array( [1/21 * ( 5 * a_s + 16 * a_s ),0,0,] ),np.array( \
[0,1/21 * ( 5 * a_s + 16 * a_s ),0,] ),np.array( [0,0,1/21 * ( 11 * \
a_s + 10 * a_s ),] ),] ),np.array( [np.array( [0,0,complex( \
0,-1/7 ) * ( 6 )**( -1/2 ) * ( a_s + -1 * a_s ),] ),np.array( \
[0,0,1/7 * ( 6 )**( -1/2 ) * ( -1 * a_s + a_s ),] ),np.array( \
[complex( 0,-1/7 ) * ( 6 )**( -1/2 ) * ( a_s + -1 * a_s ),1/7 * ( 6 \
)**( -1/2 ) * ( -1 * a_s + a_s ),0,] ),] ),np.array( [np.array( \
[1/7 * ( 2/3 )**( 1/2 ) * ( a_s + -1 * a_s ),complex( 0,-1/7 ) * ( 2/3 \
)**( 1/2 ) * ( a_s + -1 * a_s ),0,] ),np.array( [complex( 0,-1/7 ) * \
( 2/3 )**( 1/2 ) * ( a_s + -1 * a_s ),1/7 * ( 2/3 )**( 1/2 ) * ( -1 * \
a_s + a_s ),0,] ),np.array( [0,0,0,] ),] ),] ),np.array( \
[np.array( [np.array( [0,0,complex( 0,1 ) * ( 30 )**( -1/2 ) * \
( a_s + -1 * a_s ),] ),np.array( [0,0,( 30 )**( -1/2 ) * ( -1 * a_s + \
a_s ),] ),np.array( [complex( 0,1 ) * ( 30 )**( -1/2 ) * ( a_s + -1 \
* a_s ),( 30 )**( -1/2 ) * ( -1 * a_s + a_s ),0,] ),] ),np.array( \
[np.array( [0,0,0,] ),np.array( [0,0,0,] ),np.array( \
[0,0,0,] ),] ),np.array( [np.array( [0,0,0,] ),np.array( \
[0,0,0,] ),np.array( [0,0,0,] ),] ),np.array( [np.array( \
[0,0,0,] ),np.array( [0,0,0,] ),np.array( [0,0,0,] ),] \
),np.array( [np.array( [0,0,0,] ),np.array( [0,0,0,] \
),np.array( [0,0,0,] ),] ),np.array( [np.array( [1/7 * ( a_s \
+ -1 * a_s ),complex( 0,1/7 ) * ( a_s + -1 * a_s ),0,] ),np.array( \
[complex( 0,1/7 ) * ( a_s + -1 * a_s ),1/7 * ( -1 * a_s + a_s ),0,] \
),np.array( [0,0,0,] ),] ),np.array( [np.array( \
[0,0,complex( 0,1/7 ) * ( 6 )**( -1/2 ) * ( a_s + -1 * a_s ),] \
),np.array( [0,0,1/7 * ( 6 )**( -1/2 ) * ( -1 * a_s + a_s ),] \
),np.array( [complex( 0,1/7 ) * ( 6 )**( -1/2 ) * ( a_s + -1 * a_s \
),1/7 * ( 6 )**( -1/2 ) * ( -1 * a_s + a_s ),0,] ),] ),np.array( \
[np.array( [1/7 * ( 2 * a_s + 5 * a_s ),0,0,] ),np.array( [0,1/7 \
* ( 2 * a_s + 5 * a_s ),0,] ),np.array( [0,0,1/7 * ( 3 * a_s + 4 * a_s \
),] ),] ),np.array( [np.array( [0,0,complex( 0,-1/7 ) * ( a_s + \
-1 * a_s ),] ),np.array( [0,0,1/7 * ( -1 * a_s + a_s ),] \
),np.array( [complex( 0,-1/7 ) * ( a_s + -1 * a_s ),1/7 * ( -1 * a_s \
+ a_s ),0,] ),] ),] ),np.array( [np.array( [np.array( [( 30 \
)**( -1/2 ) * ( -1 * a_s + a_s ),complex( 0,-1 ) * ( 30 )**( -1/2 ) * ( \
a_s + -1 * a_s ),0,] ),np.array( [complex( 0,-1 ) * ( 30 )**( -1/2 ) \
* ( a_s + -1 * a_s ),( 30 )**( -1/2 ) * ( a_s + -1 * a_s ),0,] \
),np.array( [0,0,0,] ),] ),np.array( [np.array( [0,0,0,] \
),np.array( [0,0,0,] ),np.array( [0,0,0,] ),] ),np.array( \
[np.array( [0,0,0,] ),np.array( [0,0,0,] ),np.array( \
[0,0,0,] ),] ),np.array( [np.array( [0,0,0,] ),np.array( \
[0,0,0,] ),np.array( [0,0,0,] ),] ),np.array( [np.array( \
[0,0,0,] ),np.array( [0,0,0,] ),np.array( [0,0,0,] ),] \
),np.array( [np.array( [0,0,0,] ),np.array( [0,0,0,] \
),np.array( [0,0,0,] ),] ),np.array( [np.array( [1/7 * ( 2/3 \
)**( 1/2 ) * ( a_s + -1 * a_s ),complex( 0,1/7 ) * ( 2/3 )**( 1/2 ) * ( \
a_s + -1 * a_s ),0,] ),np.array( [complex( 0,1/7 ) * ( 2/3 )**( 1/2 \
) * ( a_s + -1 * a_s ),1/7 * ( 2/3 )**( 1/2 ) * ( -1 * a_s + a_s ),0,] \
),np.array( [0,0,0,] ),] ),np.array( [np.array( \
[0,0,complex( 0,1/7 ) * ( a_s + -1 * a_s ),] ),np.array( [0,0,1/7 * \
( -1 * a_s + a_s ),] ),np.array( [complex( 0,1/7 ) * ( a_s + -1 * a_s \
),1/7 * ( -1 * a_s + a_s ),0,] ),] ),np.array( [np.array( [1/7 * \
( 3 * a_s + 4 * a_s ),0,0,] ),np.array( [0,1/7 * ( 3 * a_s + 4 * a_s \
),0,] ),np.array( [0,0,1/7 * ( a_s + 6 * a_s ),] ),] ),] ),] )
            
            
            
#             a_23 = (2*a_p+3*a_s)
#             a_32 = (3*a_p+2*a_s)        
#             a_14 = (a_p+4*a_s)
#             a_11 = (a_p-a_s)
#             a_34 = (3*a_p+4*a_s)
#             a_16 = a_p+6*a_s
#             a_25 = 2*a_p+5*a_s
            
#             MM = np.zeros((9,9,3,3),dtype = 'complex_')
#             MM[0,0] = np.array([[a_00,0,0],[0,a_00,0],[0,0,a_00]])
#             MM[0,4] = np.array([[-a_11/np.sqrt(30),1j*a_11/np.sqrt(30),0],[-1j*a_11/np.sqrt(30),a_11/np.sqrt(30),0],[0,0,0]])
#             MM[0,5] = np.array([[0,0,-1j*a_11/np.sqrt(30)],[0,0,a_11/np.sqrt(30)],[-1j*a_11/np.sqrt(30),a_11/np.sqrt(30),0]])
#             MM[0,6] = np.array([[-a_11/(3*np.sqrt(5)),0,0],[0,-a_11/(3*np.sqrt(5)),0],[0,0,2*a_11/(3*np.sqrt(5))]])
#             MM[0,7] = np.array([[0,0,-1j*a_11/np.sqrt(30)],[0,0,-a_11/np.sqrt(30)],[-1j*a_11/np.sqrt(30),-a_11/np.sqrt(30),0]])
#             MM[0,8] = np.array([[-a_11/np.sqrt(30),1j*a_11/np.sqrt(30),0],[1j*a_11/np.sqrt(30),a_11/np.sqrt(30),0],[0,0,0]])

#             MM[1,1] = np.array([[a_23/5,0,0],[0,a_23/5,0],[0,0,a_14/5]])
#             MM[1,2] = np.array([[0,0,1j*a_11/np.sqrt(50)],[0,0,a_11/np.sqrt(50)],[1j*a_11/np.sqrt(50),a_11/np.sqrt(50),0]])
#             MM[1,3] = np.array([[a_11/5,-1j*a_11/5,0],[-1j*a_11/5,-a_11/5,0],[0,0,0]])        


#             MM[2,1] = np.array([[0,0,-1j*a_11/np.sqrt(50)],[0,0,a_11/np.sqrt(50)],[-1j*a_11/np.sqrt(50),a_11/np.sqrt(50),0]])
#             MM[2,2] = np.array([[a_14/5,0,0],[0,a_14/5,0],[0,0,a_32/5]])
#             MM[2,3] = np.array([[0,0,-1j*a_11/np.sqrt(50)],[0,0,-a_11/np.sqrt(50)],[-1j*a_11/np.sqrt(50),-a_11/np.sqrt(50),0]])

#             MM[3,1] = np.array([[a_11/5,1j*a_11/5,0],[1j*a_11/5,-a_11/5,0],[0,0,0]])
#             MM[3,2] = np.array([[0,0,1j*a_11/np.sqrt(50)],[0,0,-a_11/np.sqrt(50)],[1j*a_11/np.sqrt(50),-a_11/np.sqrt(50),0]])
#             MM[3,3] = np.array([[a_23/5,0,0],[0,a_23/5,0],[0,0,a_14/5]])   

#             MM[4,0] = np.array([[-a_11/np.sqrt(30),1j*a_11/np.sqrt(30),0],[1j*a_11/np.sqrt(30),a_11/np.sqrt(30),0],[0,0,0]])
#             MM[4,4] = np.array([[a_34/7,0,0],[0,a_34/7,0],[0,0,a_16/7]])
#             MM[4,5] = np.array([[0,0,1j*a_11/7],[0,0,a_11/7],[1j*a_11/7,a_11/7,0]])
#             MM[4,6] = np.array([[1/7*np.sqrt(2/3)*a_11,-1j*1/7*np.sqrt(2/3)*a_11,0],[-1j*1/7*np.sqrt(2/3)*a_11,-1/7*np.sqrt(2/3)*a_11,0],[0,0,0]])

#             MM[5,0] = np.array([[0,0,1j*a_11/np.sqrt(30)],[0,0,a_11/np.sqrt(30)],[1j*a_11/np.sqrt(30),a_11/np.sqrt(30),0]])
#             MM[5,4] = np.array([[0,0,-1j*a_11/7],[0,0,a_11/7],[-1j*a_11/7,a_11/7,0]])
#             MM[5,5] = np.array([[a_25/7,0,0],[0,a_25/7,0],[0,0,a_34/7]])
#             MM[5,6] = np.array([[0,0,1j*a_11/(7*np.sqrt(6))],[0,0,a_11/(7*np.sqrt(6))],[1j*a_11/(7*np.sqrt(6)),a_11/(7*np.sqrt(6)),0]])
#             MM[5,7] = np.array([[a_11/7,-1j*a_11/7,0],[-1j*a_11/7,-a_11/7,0],[0,0,0]])

#             MM[6,0] = np.array([[-a_11/(3*np.sqrt(5)),0,0],[0,-a_11/(3*np.sqrt(5)),0],[0,0,2*a_11/(3*np.sqrt(5))]])
#             MM[6,4] = np.array([[1/7*np.sqrt(2/3)*a_11,1j*1/7*np.sqrt(2/3)*a_11,0],[1j*1/7*np.sqrt(2/3)*a_11,-1/7*np.sqrt(2/3)*a_11,0],[0,0,0]])
#             MM[6,5] = np.array([[0,0,-1j*a_11/(7*np.sqrt(6))],[0,0,a_11/(7*np.sqrt(6))],[-1j*a_11/(7*np.sqrt(6)),a_11/(7*np.sqrt(6)),0]])
#             MM[6,6] = np.array([[1/21*(5*a_p+16*a_s),0,0],[0,1/21*(5*a_p+16*a_s),0],[0,0,1/21*(11*a_p+10*a_s)]])
#             MM[6,7] = np.array([[0,0,-1j*a_11/(7*np.sqrt(6))],[0,0,-a_11/(7*np.sqrt(6))],[-1j*a_11/(7*np.sqrt(6)),-a_11/(7*np.sqrt(6)),0]])
#             MM[6,8] = np.array([[1/7*np.sqrt(2/3)*a_11,-1j*1/7*np.sqrt(2/3)*a_11,0],[-1j*1/7*np.sqrt(2/3)*a_11,-1/7*np.sqrt(2/3)*a_11,0],[0,0,0]])

#             MM[7,0] = np.array([[0,0,1j*a_11/np.sqrt(30)],[0,0,-a_11/np.sqrt(30)],[1j*a_11/np.sqrt(30),-a_11/np.sqrt(30),0]])
#             MM[7,5] = np.array([[a_11/7,1j*a_11/7,0],[1j*a_11/7,-a_11/7,0],[0,0,0]])
#             MM[7,6] = np.array([[0,0,1j*a_11/(7*np.sqrt(6))],[0,0,-a_11/(7*np.sqrt(6))],[1j*a_11/(7*np.sqrt(6)),-a_11/(7*np.sqrt(6)),0]])
#             MM[7,7] = np.array([[a_25/7,0,0],[0,a_25/7,0],[0,0,a_34/7]])
#             MM[7,8] = np.array([[0,0,-1j*a_11/7],[0,0,-a_11/7],[-1j*a_11/7,-a_11/7,0]])

#             MM[8,0] = np.array([[-a_11/np.sqrt(30),-1j*a_11/np.sqrt(30),0],[-1j*a_11/np.sqrt(30),a_11/np.sqrt(30),0],[0,0,0]])
#             MM[8,6] = np.array([[1/7*np.sqrt(2/3)*a_11,1j*1/7*np.sqrt(2/3)*a_11,0],[1j*1/7*np.sqrt(2/3)*a_11,-1/7*np.sqrt(2/3)*a_11,0],[0,0,0]])
#             MM[8,7] = np.array([[0,0,1j*a_11/7],[0,0,-a_11/7],[1j*a_11/7,-a_11/7,0]])
#             MM[8,8] = np.array([[a_34/7,0,0],[0,a_34/7,0],[0,0,a_16/7]])
        else:
            print('mising AC stark term!')
        
        MM = np.einsum('ijkl->klij',MM)    
        if pol.shape == (3,1):
            P = np.einsum('ijkl,jm->ikl',MM,pol)
            Q = np.einsum('ij,jlm->ilm',np.conjugate(pol.T),P)
        else:
            P = np.einsum('ijkl,jm->ikl',MM,pol.T)
            Q = np.einsum('ij,jlm->ilm',np.conjugate(pol),P)
        Q = Q / np.linalg.norm(pol)
        mat = rot_state_op(b,Q[0,:],a_00)
        
        return delta_spec*mat*E00/a_00     

    H_ls_uc = -H_ls(b_uc,p['alpha_p'],p['alpha_s'],1,p['trap_pol'])

    if np.linalg.norm(p['E_dc']) > 0.001:
        M_r2,M_z,M_r1 = dm_ops_sph(b_uc)
        M_x = (M_r1 + np.conjugate(M_r2.T))/np.sqrt(2)
        M_y = (M_r1 - np.conjugate(M_r2.T))*1j/np.sqrt(2)
        if len(p['E_dc'].shape) == 1:
            H_S_uc = - p['muE'] * (0 * M_x + 0 * M_y + p['E_dc'][:, np.newaxis, np.newaxis] * M_z)    
        else:
            Mxyz = np.array([M_x,M_y,M_z])
            H_S_uc = - p['muE'] * np.inner(p['E_dc'].transpose(),np.einsum('ijk->jki',Mxyz))
#             H_S_uc = - p['muE'] * (p['E_dc'][0] * M_x + p['E_dc'][1] * M_y + p['E_dc'][2] * M_z)
#         H_S_uc = - p['muE'] * (0 * M_x + 0 * M_y + p['E_dc'][:, np.newaxis, np.newaxis] * M_z)    
    else:
        H_S_uc = 0
    
    if not(p['B'].shape[0] > 1 & p['E_td_00'].shape[0] > 1):
        H_B_tot_uc = (H_rot_uc + H_hf_sr_uc + H_hf_ss_uc + H_hf_eQ_uc + H_ls_uc*p['E_td_00'][:, np.newaxis, np.newaxis] +
            (H_z_rot_uc + H_z_nuc_uc) * p['B'][:, np.newaxis, np.newaxis]) + H_S_uc
    else:
        raise ValueError('cannot scan multiple axes')
        
    return np.real(H_B_tot_uc), b_uc


def map2basis(Tv):
    # finds the "closest" uncoupled basis index value
    M = np.zeros(Tv.shape)
    inds = []
    for j in range(Tv.shape[0]):
        v = Tv[:,j]
        v_amp = np.abs(v)
        ind = np.argmax(v_amp)
        M[ind,j] = 1
        inds.append(ind)
    return M, np.array(inds)

def disp_vec(b_uc, Evv, qn_uc, b_print=True):
    # display components corresponding to uc basis
    if type(qn_uc) == list:
        _,inds_vind2uc = map2basis(Evv)
        uc_ind = b_uc.get_ind(*qn_uc)[0]
        ev_ind = np.where(inds_vind2uc==uc_ind)[0][0]
        vec_amps = Evv[:,ev_ind] * Evv[:,ev_ind]
        vec_amps0 = Evv[:,ev_ind]
        amp_argsort = np.squeeze(np.argsort(-vec_amps,axis=0))
        pos = np.argmax(np.squeeze(-np.sort(-vec_amps,axis=0)) < 0.0001)
        ddf = b_uc.v.loc[amp_argsort[0:pos],['m_i1','m_i2','n','m_n']]
        ddf['amps'] =  np.squeeze(vec_amps0[amp_argsort[0:pos]])
    elif type(qn_uc) == int:
        vec_amps = Evv[:,qn_uc] * Evv[:,qn_uc]        
        vec_amps0 = Evv[:,qn_uc]
        amp_argsort = np.squeeze(np.argsort(-vec_amps,axis=0))
        pos = np.argmax(np.squeeze(-np.sort(-vec_amps,axis=0)) < 0.0001)
        ddf = b_uc.v.loc[amp_argsort[0:pos],['m_i1','m_i2','n','m_n']]
        ddf['amps'] =  np.squeeze(vec_amps0[amp_argsort[0:pos]])       
        ev_ind = qn_uc
    if b_print:
        print('Ev index: ' + str(ev_ind))    
        print(ddf)
    return ev_ind

def get_max_ind(v):
    return np.argmax(np.abs(v),axis=0)

def calc_eig(H):
    E0, Ev0 = np.linalg.eigh(H)
    E0 = np.real(E0)
    
    if Ev0.shape[0] == 1:
        Ev0 = np.squeeze(Ev0)
    else:
        sort_ind = np.argsort(E0,axis=1)
        E0 = np.take_along_axis(E0,sort_ind,axis=1)
        Ev0 = np.take_along_axis(Ev0,sort_ind[:,None,:],axis=2)
    return E0, Ev0

def calc_uwave(E,Ev,b_uc,ss):
    
    def map2roteig(Tv, b_uc):
        inds = np.zeros(Tv.shape[0])
        for j in range(Tv.shape[0]):
            v = Tv[:,j]
            v_amp = np.abs(v)
            ind = np.argmax(v_amp)
            if b_uc.v.loc[ind,'m_n'] == 0:
                inds[j] = 0
            elif b_uc.v.loc[ind,'m_n'] == 1 or b_uc.v.loc[ind,'m_n'] == -1:
                v_amp[ind] = 0
                ind2 = np.argmax(v_amp)
                if np.sign(v[ind])*np.sign(v[ind2]) > 0:
                    inds[j] = 1
                else:
                    inds[j] = -1
        return inds
    
    sind = disp_vec(b_uc, Ev, ss)    
    
    H_dm_uc = dm_ops_sph(b_uc)
    H_dm_v = []
    for i in range(3):
        H_dm_v.append(np.transpose(Ev)@H_dm_uc[i]@Ev)
                
    thresh = 0.00001
    all_ind = np.array([],dtype=np.int8)
    all_ind_uc = np.array([],dtype=np.int8)
    vals = np.empty([3],dtype=np.int8)
    vals_norm = np.array([],dtype=np.int8)
    pols = ['sig-','pi','sig+']
    trans_type = []

    mat_el = (H_dm_v[0][sind,:]**2+H_dm_v[1][sind,:]**2+H_dm_v[2][sind,:]**2)
    ind = (np.argwhere(np.abs(mat_el) > thresh))
    all_ind = np.append(all_ind,ind)
    all_ind_uc = np.append(all_ind_uc,get_max_ind(Ev[:,ind]))    
    mat_el_all = np.array([H_dm_v[0][sind,:],H_dm_v[1][sind,:],H_dm_v[2][sind,:]])
    for i in ind:
        vals = np.vstack([vals, np.squeeze(mat_el_all[:,i])])
        vals_norm = np.append(vals_norm,mat_el[i])
    vals = vals[1:,:]
    
    rot_inds = map2roteig(Ev, b_uc)

    df = b_uc.v.loc[all_ind_uc,['m_i1','m_i2','n','m_n']]
    # df['m_n'] = df['m_n']
    df['mF'] = df[['m_i1','m_i2','m_n']].sum(axis=1)
#     df['rot'] = pd.Series(rot_inds)
    df[['sig-','pi','sig+']] = pd.DataFrame(index=df.index,data=vals).applymap('{:,.5f}'.format)
    df['mat_el_abs'] = np.abs(vals_norm)
    # pols = ['sig-','pi','sig+']
    #df['pol'] = trans_type

    
    df['E [kHz]'] = E[0,all_ind] - E[0,all_ind[0]]
    df['rot'] = rot_inds[all_ind]
    df['v_ind'] = all_ind
    print(df.sort_values(['mat_el_abs'],ascending=False).drop(['mat_el_abs'],axis=1))
    #df['cost'] = df['cost'].map('${:,.2f}'.format)
    #print(df)

    
def calc_uwave_vals(E,Ev,b_uc,ss):
    
    def map2roteig(Tv, b_uc):
        inds = np.zeros(Tv.shape[0])
        for j in range(Tv.shape[0]):
            v = Tv[:,j]
            v_amp = np.abs(v)
            ind = np.argmax(v_amp)
            if b_uc.v.loc[ind,'m_n'] == 0:
                inds[j] = 0
            elif b_uc.v.loc[ind,'m_n'] == 1 or b_uc.v.loc[ind,'m_n'] == -1:
                v_amp[ind] = 0
                ind2 = np.argmax(v_amp)
                if np.sign(v[ind])*np.sign(v[ind2]) > 0:
                    inds[j] = 1
                else:
                    inds[j] = -1
        return inds
    
    sind = disp_vec(b_uc, Ev[0], ss, b_print=False)    
    
    H_dm_uc = dm_ops_sph(b_uc)
    H_dm_v_all = []
    for i in range(3):
        H_dm_v = np.zeros(Ev.shape[0:2])
        for n, vv in enumerate(Ev):
            H_dm_v[n,:] = (np.transpose(vv)@H_dm_uc[i]@vv)[sind,:]
        H_dm_v_all.append(H_dm_v)

    return H_dm_v_all # H_dm_v[:][sind,:]

def calc_uwave_abs(E,Ev,b_uc,ss):
    
    def map2roteig(Tv, b_uc):
        inds = np.zeros(Tv.shape[0])
        for j in range(Tv.shape[0]):
            v = Tv[:,j]
            v_amp = np.abs(v)
            ind = np.argmax(v_amp)
            if b_uc.v.loc[ind,'m_n'] == 0:
                inds[j] = 0
            elif b_uc.v.loc[ind,'m_n'] == 1 or b_uc.v.loc[ind,'m_n'] == -1:
                v_amp[ind] = 0
                ind2 = np.argmax(v_amp)
                if np.sign(v[ind])*np.sign(v[ind2]) > 0:
                    inds[j] = 1
                else:
                    inds[j] = -1
        return inds
    
    sind = disp_vec(b_uc, Ev, ss)    
    
    H_dm_uc = dm_ops_sph(b_uc)
    H_dm_v = []
    for i in range(3):
        H_dm_v.append(np.transpose(Ev)@H_dm_uc[i]@Ev)
                
    thresh = 0.00001
    all_ind = np.array([],dtype=np.int8)
    all_ind_uc = np.array([],dtype=np.int8)
    vals = np.empty([3],dtype=np.int8)
    vals_norm = np.array([],dtype=np.int8)
    pols = ['sig-','pi','sig+']
    trans_type = []

    mat_el = (H_dm_v[0][sind,:]**2+H_dm_v[1][sind,:]**2+H_dm_v[2][sind,:]**2)
    ind = (np.argwhere(np.abs(mat_el) > thresh))
    all_ind = np.append(all_ind,ind)
    all_ind_uc = np.append(all_ind_uc,get_max_ind(Ev[:,ind]))    
    mat_el_all = np.array([H_dm_v[0][sind,:],H_dm_v[1][sind,:],H_dm_v[2][sind,:]])
    for i in ind:
        vals = np.vstack([vals, np.squeeze(mat_el_all[:,i])])
        vals_norm = np.append(vals_norm,mat_el[i])
    vals = vals[1:,:]
    
    rot_inds = map2roteig(Ev, b_uc)

    df = b_uc.v.loc[all_ind_uc,['m_i1','m_i2','n','m_n']]
    # df['m_n'] = df['m_n']
    df['mF'] = df[['m_i1','m_i2','m_n']].sum(axis=1)
#     df['rot'] = pd.Series(rot_inds)
    df[['sig-','pi','sig+']] = pd.DataFrame(index=df.index,data=vals).applymap('{:,.5f}'.format)
    df['mat_el_abs'] = np.abs(vals_norm)
    # pols = ['sig-','pi','sig+']
    #df['pol'] = trans_type

    
    df['E [kHz]'] = E[0,all_ind] - 2*p0['Bv'] - E[0,1]
    df['rot'] = rot_inds[all_ind]
    df['v_ind'] = all_ind
    print(df.sort_values(['mat_el_abs'],ascending=False).drop(['mat_el_abs'],axis=1))
    #df['cost'] = df['cost'].map('${:,.2f}'.format)
    #print(df)
    
    
def color_map_blue(norm_val):
    return (1-np.exp(np.log(norm_val)),1-np.exp(np.log(norm_val)),1)
def color_map_red(norm_val):
    return (1,1-np.power(norm_val,1),1-np.power(norm_val,1))
