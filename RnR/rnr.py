# rewriting Ryan's RNR script (in the end, i didnt change much)
# question: what about gravity effect?
# following Lee's thesis, page 187

import numpy as np
import scipy.constants as const
import time
from scipy.linalg import expm
# import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from scipy.stats import rv_discrete
from scipy.stats import uniform

# constants
c = const.c
e_0 = const.epsilon_0
h = const.h
hbar = const.hbar
kB = const.k
a_0 = const.physical_constants["atomic unit of length"][0]
amu = const.physical_constants["atomic mass constant"][0]
m_cs = 133*amu
m_nacs = 155.8952*amu
kb = 1.381*(10**-23) #Boltzmann constant in J/K

def get_nbar(f_trap,T):
    # f_trap: trapping freq in Hz
    omega_trap = 2*np.pi*f_trap
    boltz_factor = np.exp(hbar*omega_trap/(kB*T))
    return 1/(boltz_factor-1)

def Pn(n,nbar):
    # probability of n given nbar
    return (nbar**n)/((1+nbar)**(n+1))

def step(n, trap_freq, m, theta, t_rr):
    Etotal = (n+0.5)*hbar*(2*np.pi*trap_freq)
    # KE = Etotal*(np.sin(theta)**2)
    # PE = Etotal*(np.sin(cos)**2)
    x_init = np.sqrt(2*Etotal/m)*np.sin(theta)/(2*np.pi*trap_freq)
    v_init = np.sqrt(2*Etotal/m)*np.cos(theta)
    x_fin = x_init + v_init*t_rr

    return np.array([v_init,x_fin])

class thermal_quantum_gen(rv_discrete):
    def _pmf(self, n, n_bar):
        return Pn(n,n_bar) #eqn c.7.2

thermal_quantum = thermal_quantum_gen(name="thermal_quantum")

def waist(w0, lbda, z):
    w = w0*np.sqrt(1+((lbda*z)/(np.pi*(w0**2)))**2)
    return w

def montecarlo(trap_freq_arr, nb_arr, m, t_rr, trap_depth, lbda_tw, reps):

    # t_rr release and recapture time

    p_survive = 0
    #nb_arr = [get_nbar(trap_freq_arr[0], temps_arr[0]),get_nbar(trap_freq_arr[1], temps_arr[1]),get_nbar(trap_freq_arr[2], temps_arr[2])]
    #nb_arr = [.412, .103, .049]
    velocity_arr_ax = []
    velocity_arr_r1 = []
    velocity_arr_r2 = []
    energy_arr = []
    
    # Equations for beam waist
    #w0_tw = np.sqrt(kb*trap_depth/((trap_freq_arr[1]**2)*m*(np.pi**2)))
    # use the ratio between axial and radial instead
    w0_tw = lbda_tw*(trap_freq_arr[1]/trap_freq_arr[0])/(np.pi*np.sqrt(2)) 
    
    for i in range(0,reps):
        n_arr = [thermal_quantum.rvs(n_bar=nb_arr[0]),thermal_quantum.rvs(n_bar=nb_arr[1]),thermal_quantum.rvs(n_bar=nb_arr[2])]
        theta_arr = [uniform(scale=2*np.pi).rvs(),uniform(scale=2*np.pi).rvs(),uniform(scale=2*np.pi).rvs()]
        final_conditions_arr = [step(n_arr[0], trap_freq_arr[0], m, theta_arr[0], t_rr),step(n_arr[1], trap_freq_arr[1], m, theta_arr[1], t_rr),step(n_arr[2], trap_freq_arr[2], m, theta_arr[2], t_rr)]

        # Compute final kinetic and potential energies
        k_final = np.array([(1/2)*m*(final_conditions_arr[0][0])**2, (1/2)*m*(final_conditions_arr[1][0])**2, (1/2)*m*(final_conditions_arr[2][0])**2])

        # u_final = kb*trap_depth*(1-(w0_tw**2)/(waist(w0_tw, lbda_tw, final_conditions_arr[0][1])**2)*np.exp(-2*(final_conditions_arr[1][1]**2 + final_conditions_arr[2][1]**2)/(waist(w0_tw, lbda_tw, final_conditions_arr[0][1])**2)))        

        u_final = kb*trap_depth*(1-(w0_tw**2)/(waist(w0_tw, lbda_tw, final_conditions_arr[0][1])**2)*np.exp(-2*(final_conditions_arr[1][1]**2 + final_conditions_arr[2][1]**2)/(waist(w0_tw, lbda_tw, final_conditions_arr[0][1])**2)))    


        # print('final energy arr is ', final_energy_arr)
        # print('trap depth energy is ', kb*trap_depth)
        e_final = np.sum(k_final) + u_final
                                 
        #velocity_arr_ax.append(abs(final_energy_arr[0][2]))
        #velocity_arr_r1.append(abs(final_energy_arr[1][2]))
        #velocity_arr_r2.append(abs(final_energy_arr[2][2]))
        energy_arr.append(e_final)
        
        if (e_final <= kb*trap_depth):
            #p_survive += 1
            #Account for detection mismatches. In this case it's a detection of ground as loss
            if(uniform(scale=1).rvs()>.04):
                p_survive += 1
        #else:
            #Now it's a loss detected as ground
            #if(uniform(scale=1).rvs()<.13):
                #p_survive += 1
    p_survive = p_survive/reps
    #plt.hist(velocity_arr_ax, 40)
    #plt.title('Axial velocity distribution at t=%.7f s' % t_rr)
    #plt.show()
    #plt.hist(velocity_arr_r1, 40)
    #plt.title('Radial 1 velocity distribution at t=%.7f s' % t_rr)
    #plt.show()
    #plt.hist(velocity_arr_r2, 40)
    #plt.title('Radial 2 velocity distribution at t=%.7f s' % t_rr)
    #plt.show()
    #plt.hist(energy_arr, 40)
    #plt.title('Final energy distribution at t=%.7f s' % t_rr)
    #plt.vlines(kb*trap_depth, 0, 40)
    #plt.show()
    return p_survive


if __name__ == '__main__':

    thermal_quantum = thermal_quantum_gen(name="thermal_quantum")

    # temperature = 150e-6
    # f_trap = 150e3

    # nbar = get_nbar(f_trap,temperature)


    # n = np.linspace(0,200,201)
    # # n2 = np.linspace(1,100,1000)


    # # print(np.sum(Pn(n,nbar))*(n[2]-n[1])+np.sum(Pn(n2,nbar))*(n2[2]-n2[1]))
    # print(np.sum(Pn(n,nbar)))

    # # sample from thermal quantum
    # samples = np.zeros(10000)
    # for i in range(10000):
    #     samples[i] = thermal_quantum.rvs(nbar)

    ### PLOT TO TEST THE DISTRUBTION
    # plt.figure()
    # counts, bins = np.histogram(samples,bins=5000)
    # print("numb of counts",np.sum(counts))
    # print("summing prob dis:",np.sum(counts)/10000)
    # plt.stairs(counts/10000,bins)
    # # plt.show()

    # plt.figure()
    # plt.plot(n,Pn(n,nbar),'.')
    # # plt.plot(n2,Pn(n2,nbar))


    #data from RSC release_recapture
    survival_noRearr_rr = np.array([1.0000,    0.9451,    0.9050,    0.8183,    0.9248,    0.7090,    0.4331,    0.2234,    0.1667])
    survival_noRearr_rr_err = np.array([0.0970,    0.0916,    0.0921,    0.0882,    0.0918,    0.0816,    0.0674,    0.0471,    0.0415])
    survival_Rearr_rr = np.array([1.0000,    0.8393,    1.0282,    1.0333,    0.7286,    0.6602,    0.3243,    0.1649,    0.1636])
    survival_Rearr_rr_err = np.array([0.0815,    0.0764,    0.0827,    0.0846,    0.0727,    0.0684,    0.0492,    0.0350,    0.0347])
    params_rr = 1e-6*np.array([1,   4,    7,    10,    50,    100,    200,    400,    500])
    
    trap_depth_Hz = 1.3867e6
    trap_depth = trap_depth_Hz*h/kB #Convert from trap depth in Hz to mK
    trap_waist = 1.19e-6 #Use our previously calibrated value
    lambda_trap = 1064e-9
    rFreq = np.sqrt(4*trap_depth_Hz*h/m_nacs/trap_waist**2)/2/np.pi
    axFreq = np.sqrt(2*trap_depth_Hz*h/m_nacs/(np.pi*trap_waist**2/lambda_trap)**2)/2/np.pi
    trap_freqs = [axFreq, rFreq,rFreq] 

    reps = 200
    #times = np.linspace(0,500*(10**-6),20)
    times  = params_rr
    temps = np.linspace(2e-9,200e-9,5)

    nbRad = 0.15
    #nbAx = 0.15
    nbArr = np.linspace(0.1,6,10)

    mse_Rearr = []
    mse_noRearr = []

    i=0
    for nbAx in nbArr:
        print('Testing axial nb ', nbAx)
        outcomes = []
        for t_rr in times:
            outcomes.append(montecarlo(trap_freqs, [nbAx,nbRad,nbRad] , m_nacs, t_rr, trap_depth, lambda_trap, reps)) 
        plt.plot(1e6*times, outcomes, label=['Axial nBar = ',nbAx], linestyle='--')
        mse_Rearr.append(np.sum((np.array(outcomes) - survival_Rearr_rr)**2))
        mse_noRearr.append(np.sum((np.array(outcomes) - survival_noRearr_rr)**2))
        i += 1
    # for nbRad in nbArr:
    #     print('Testing radial nb ', nbRad)
    #     outcomes = []
    #     for t_rr in times:
    #         outcomes.append(montecarlo(trap_freqs, [nbAx,nbRad,nbRad] , m_nacs, t_rr, trap_depth, lambda_trap, reps)) 
    #     plt.plot(1e6*times, outcomes, label=['Radial nBar = ',nbRad], linestyle='--')
    #     mse_Rearr.append(np.sum((np.array(outcomes) - survival_Rearr_rr)**2))
    #     mse_noRearr.append(np.sum((np.array(outcomes) - survival_noRearr_rr)**2))
    #     i += 1
    plt.errorbar(1e6*params_rr, survival_noRearr_rr,survival_noRearr_rr_err, label='data no rearr')
    plt.errorbar(1e6*params_rr, survival_Rearr_rr,survival_Rearr_rr_err, label='data rearr')
    plt.xlabel('Release time (us)')
    plt.ylabel('Recapture probability')
    leg = plt.legend()

    plt.figure()
    plt.plot(nbArr,mse_noRearr, label='no Rearr')
    plt.plot(nbArr,mse_Rearr, label='Rearr')
    plt.xlabel('nBar')
    plt.ylabel('MSE')
    plt.legend()   

    plt.show()

