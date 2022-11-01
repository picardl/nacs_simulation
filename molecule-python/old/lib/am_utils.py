import numpy as np
import pandas as pd
import py3nj

def arange(val1, val2=None, ax=0):
    if val2 is None:
        val2 = -val1
    val_min = np.min([val1,val2])
    val_max = np.max([val1,val2])
    return np.arange(val_min, val_max+1, 1)

class am_basis():

    def __init__(self, cols=None, arr=None): 
        if arr is None:
            self.v = self.build_basis(cols)        
        else:
            self.v = pd.DataFrame(arr, columns=cols)
        self.set_vals()
        
    def build_basis(self, cols):    
        # cols contains list of lists
        # [['i_na', 3/2, 1], ['s_na', 1/2, 1]]        
        col_names = []
        grid_list = []    
        for c in cols:
            col_names.append(c[0])
            grid_list.append(c[1])        
            if c[2]:            
                col_names.append('m_'+c[0]) 
                gt = []
                c_f = c[1]
                if not isinstance(c[1], list):
                    c_f = [c[1]]                  
                for f in c_f:
                    gt.extend(arange(-f,f).tolist())
                grid_list.append(gt)

        grid_arr = np.meshgrid(*grid_list)

        n_col = len(col_names)
        n_row = grid_arr[0].size
        d_np = np.empty((n_row,n_col))        
        for i,g in enumerate(grid_arr):
            d_np[:,i] = g.flatten()

        df = pd.DataFrame(d_np,columns=col_names)
        for j in col_names:
            if 'm_'+j in col_names:
                df = df[(df['m_'+j]<=df[j]) & (df['m_'+j]>=-df[j])]
            df.drop_duplicates(inplace=True,ignore_index=True)
        
            
        return df # pd.DataFrame(d_np,columns=col_names)

    def set_vals(self):
        self.cols = self.v.columns.to_list()
        self.nrow = self.v.shape[0]
        self.ncol = self.v.shape[1]
        self.np = self.v.to_numpy()   
        self.v2 = (2*self.v).astype('int16')
        # self.op = basic_ops(self)
        
    def get_cols(self):
        return self.v.columns.to_list()
        
    def disp(self):
        print(self.v)
    
    def get_ind(self, *argv):
        i = 0
        for arg in argv:
            if i == 0:
                cond = (self.v[arg[0]]==arg[1])
                i = 1
            else:
                cond = cond & (self.v[arg[0]]==arg[1])
        return self.v.index[cond].values

        
def join_basis(b_list):
    col_names = []
    b_list_nrow = []
    for b in b_list:
        for c in b.cols:
            col_names.append(c) 
        b_list_nrow.append(range(b.nrow))
    
    g_ind = np.meshgrid(*b_list_nrow)
    b_concat = []
    for b,g in zip(b_list,g_ind):
        b_concat.append(b.np[g.flatten()])
    b_np = np.concatenate(b_concat, axis=1)
    return am_basis(col_names,b_np)


def b_op(f,b,is_py3nj,*arg):
    x,y = np.meshgrid(range(0,b.nrow),range(0,b.nrow))
    x = np.reshape(x, (-1))
    y = np.reshape(y, (-1))
    
    if is_py3nj:
        bx = b.v2.iloc[x].reset_index(drop=True)
        by = b.v2.iloc[y].reset_index(drop=True)
        return np.reshape(f(bx,by,*arg),(b.nrow,b.nrow))
    else:
        bx = b.v.iloc[x].reset_index(drop=True)
        by = b.v.iloc[y].reset_index(drop=True)      
        return np.reshape(f(bx,by,*arg).to_numpy(),(b.nrow,b.nrow))
    
def sph_ops(b,j):
    _,_,delta_spec = get_spectators(b,[j,'m_'+j],mat=True)
    f1 = b_op(lambda b1,b2,j: (-1)**(b1[j]-b1['m_'+j]), b, 0,j)    
    f2 = b_op(lambda b1,b2,j: np.sqrt(b2[j]*(b2[j]+1)*(2*b2[j]+1)),b,0,j)
    delta_j = b_op(lambda b1,b2,j: b1[j]==b2[j],b,0,j)
    cg_z = b_op(lambda b1,b2,j,q: py3nj.wigner3j(b1[j],    int(2*1),   b2[j],
                                              -b1['m_'+j], int(2*q), b2['m_'+j]), b,1,j,0)
    cg_r1 = b_op(lambda b1,b2,j,q: py3nj.wigner3j(b1[j],   int(2*1),   b2[j],
                                              -b1['m_'+j], int(2*q), b2['m_'+j]), b,1,j,1)
    cg_r2 = b_op(lambda b1,b2,j,q: py3nj.wigner3j(b1[j],       int(2*1), b2[j],
                                              -b1['m_'+j], int(2*q), b2['m_'+j]), b,1,j,-1)
    cg_y = 1j/2*np.sqrt(2) * (cg_r1+cg_r2)
    cg_x = 1/2*np.sqrt(2) * (cg_r2-cg_r1)

    M_z = delta_spec*f1*f2*delta_j*cg_z
    M_x = delta_spec*f1*f2*delta_j*cg_x
    M_y = delta_spec*f1*f2*delta_j*cg_y
    return (M_x,M_y,M_z)
    
def basic_ops(b):
    # find pairs
    op_dict = {};    
    cols = b.get_cols()
    
    # spherical terms
    j_list = []
    for c in cols:
        if 'm_'+c in cols:
            j_list.append(c)    
    for j in j_list:        
        jx,jy,jz = sph_ops(b,j)
        op_dict[j+'_x'] = jx
        op_dict[j+'_y'] = jy
        op_dict[j+'_z'] = jz
        op_dict[j+'^2'] = np.real(jx@jx + jy@jy + jz@jz)

    # delta funcs
    for j in cols:
        op_dict['eq_'+j] = b_op(lambda b1,b2,j: b1[j]==b2[j],b,0,j)
        
    return op_dict

# parity operation
def par_op(f, b):
    return b_op(lambda b1,b2: (-1)**(f(b1,b2)), b, 0)        

def scal_op(f, b):
    return b_op(lambda b1,b2: f(b1,b2),b,0)

def dot_op(b, j1, j2):
    return np.real(b.op[j1+'_x']@b.op[j2+'_x']+b.op[j1+'_y']@b.op[j2+'_y']+b.op[j1+'_z']@b.op[j2+'_z'])  
def diag_op(b, j):
    v = b.v[j]
    return np.diag(v)

def op(b,op_str):
    if op_str[-2:] == '^2':
        J = diag_op(b, op_str[0:-2])
        return J*(J+1)
    if op_str[0:2] == 'm_':
        j = op_str[2:]
        return diag_op(b, 'm_'+j)        
    ind_dot = op_str.find('*')  
    if ind_dot > -1:
        ind_com = op_str.find(',')
        J1 = diag_op(b, op_str[0:ind_dot])
        J2 = diag_op(b, op_str[(ind_dot+1):ind_com])
        J = diag_op(b, op_str[(ind_com+1):])
        return (J*(J+1) - J1*(J1+1) - J2*(J2+1))/2

def couple_basis_single(b,j1,j2,j):

    cols_ex = ['m_'+j1,'m_'+j2]
    cols_in = get_spectators(b,cols_ex)[1] + [j,'m_'+j]
    
    bv = b.v.drop(cols_ex, axis=1, inplace=False).drop_duplicates(inplace=False)
    b_concat = []
    for i,r in bv.iterrows():
        j_min = np.abs(r[j1]-r[j2])
        j_max = np.abs(r[j1]+r[j2])
        for jval in arange(j_min,j_max):
            jvals = np.atleast_2d(arange(jval)).T
            s = jvals.shape
            b_concat.append(np.concatenate([np.tile(np.array(r),s),jval*np.ones(s),jvals],axis=1))
            
    bc = am_basis(cols_in,np.concatenate(b_concat,axis=0))
    M = trans_matrix(b,bc,j1,j2,j)
    return (bc, M)

def trans_matrix(b1,b2,j1,j2,j):
    
    # transformation matrix from b1 (uncoupled) to b2 (coupled)
    x,y = np.meshgrid(range(0,b1.nrow),range(0,b1.nrow))
    x = np.reshape(x, (-1))
    y = np.reshape(y, (-1))
    
    bx = b1.v2.iloc[x].reset_index(drop=True)
    by = b2.v2.iloc[y].reset_index(drop=True)
    #cg = py3nj.wigner3j(bx[j1],bx[j2],by[j],
    #                bx['m_'+j1],bx['m_'+j2],-by['m_'+j])*np.sqrt(by[j]+1)
    cg = py3nj.clebsch_gordan(bx[j1],bx[j2],by[j],
                       bx['m_'+j1],bx['m_'+j2],by['m_'+j])

    # spectator columns
    b1_s = get_spectators(b1,['m_'+j1,'m_'+j2])[0].iloc[x].reset_index(drop=True)
    b2_s = get_spectators(b2,[j,'m_'+j])[0].iloc[y].reset_index(drop=True)
    del_spec = np.prod(np.equal(b1_s,b2_s),axis=1)
    
    return np.reshape(np.expand_dims(cg*del_spec,1),(b1.nrow,b1.nrow))                       


def couple_basis(b, j=None, j_list=None):
    if j:
        return couple_basis_single(b,*j)
    else:
        bc = b
        MM = None
        for jj in j_list:
            bc, M = couple_basis_single(bc,*jj)
            if MM is None:
                MM = M
            else:
                MM = M @ MM
        return bc, MM
        
def get_spectators(b,j_list,mat=False):
    cols = b.get_cols()
    for j in j_list:
        cols.remove(j)
    if mat:
        M = np.zeros([b.nrow,b.nrow])
        v_spec = b.v[cols].copy()
        for i,r in v_spec.iterrows():
            M[i,:] = np.prod(np.equal(v_spec, r), axis=1)  
        return (v_spec, cols, M)
    else:
        return (b.v[cols].copy(), cols)        
