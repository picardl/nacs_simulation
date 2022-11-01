import scipy as sp
import scipy.constants

hbar = sp.constants.hbar
h = sp.constants.h # m^2 kg s^-1 = J s

muB = 9.274009994e-24 # sp.constants.physical_constants["Bohr magneton"][0]  # J T^-1
muN = 5.050783699e-27 # sp.constants.physical_constants["neutron mag. mom."][0] # J T^-1

ge = sp.constants.physical_constants["electron g factor"][0]

A_na = h * 885.81306440e6 # Hz
A_na_MHz = 885.81306440

A_cs = h * 2.2981579425e9 # Hz
A_cs_MHz = 2.2981579425e3  

gS = 2.0023193043622
gI_na = -0.00080461080
gI_cs = -0.00039885395

e = sp.constants.e # 1.60217663e-19 # [C]
a0 = sp.constants.physical_constants['Bohr radius'][0]
me = sp.constants.m_e
eps0 = sp.constants.epsilon_0
Eh = sp.constants.physical_constants['Hartree energy'][0]

c = sp.constants.c    
    
