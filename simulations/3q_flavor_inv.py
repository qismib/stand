# ---------------------------- 3q_flavor_inv.py ----------------------------
#---------------------------------------------------------------------------
# Reproduce plots from https://arxiv.org/pdf/2207.03189 figure 2 for N=3.

import numpy as np
import qutip as qt
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg')
from scipy.constants import c, elementary_charge





def main():
    N = 3 # number of neutrinos, [N] = 1
    mu = 10e-4 * N / 4 / 10e6 # mu = (delta m)^2 * N/4E, where delta m: mass difference between mass eigenstates, E: typical neutrino energies (10 MeV) TODO check units!!!
    theta_nu = 0.195

    # g12 = mu * (1 - np.cos(np.arccos(0.9)*abs(1-2)/(N-1)))
    # g13 = mu * (1 - np.cos(np.arccos(0.9)*abs(1-3)/(N-1)))
    # g23 = mu * (1 - np.cos(np.arccos(0.9)*abs(2-3)/(N-1)))

    g12 = 1e-1
    g13 = 1e-1
    g23 = 1e-2

    tlist = np.linspace(0, 100, 10000)

    # one-qubit operators
    sx1 = qt.tensor(qt.sigmax(), qt.qeye(2), qt.qeye(2))
    sy1 = qt.tensor(qt.sigmay(), qt.qeye(2), qt.qeye(2))
    sz1 = qt.tensor(qt.sigmaz(), qt.qeye(2), qt.qeye(2))

    sx2 = qt.tensor(qt.qeye(2), qt.sigmax(), qt.qeye(2))
    sy2 = qt.tensor(qt.qeye(2), qt.sigmay(), qt.qeye(2))
    sz2 = qt.tensor(qt.qeye(2), qt.sigmaz(), qt.qeye(2))

    sx3 = qt.tensor(qt.qeye(2), qt.qeye(2), qt.sigmax())
    sy3 = qt.tensor(qt.qeye(2), qt.qeye(2), qt.sigmay())
    sz3 = qt.tensor(qt.qeye(2), qt.qeye(2), qt.sigmaz())

    # initial state
    psi0 = qt.tensor(qt.basis(2,0), qt.basis(2,1), qt.basis(2,0)) # |010>

    # Hamiltonians TODO add free hamiltonian to take into account all flavor changing processes
    H_Heisenberg = (np.sin(2*theta_nu)*sx1 - np.cos(2*theta_nu)*sz1) + (np.sin(2*theta_nu)*sx2 - np.cos(2*theta_nu)*sz2) + (np.sin(2*theta_nu)*sx3 - np.cos(2*theta_nu)*sz3) # free Hamiltonian
    H_Heisenberg = g12*(sx1*sx2 + sy1*sy2 + sz1*sz2) + g13*(sx1*sx3 + sy1*sy3 + sz1*sz3) + g23*(sx2*sx3 + sy2*sy3 + sz2*sz3) # interaction Hamiltonian

    H_Trotter1 = H_Heisenberg
    H_Trotter1  = 0.25 * (sx1*sz2*sy3 + sz1*sx2*sy3) * (g12*g13+g12*g23-3*g13*g23)
    H_Trotter1 += 0.25 * (sx1*sy2*sz3 + sz1*sy2*sz3) * (g12*g13+g13*g23-3*g12*g23)
    H_Trotter1 += 0.25 * (sy1*sx2*sz3 + sy1*sz2*sx3) * (g12*g23+g13*g23-3*g12*g13)

    H_Trotter2  = H_Heisenberg
    H_Trotter2 += sx1*sx2 * (13*g12*g13**2/48 - 11*g13**2*g23/48 - 1*g12*g13*g23/24 + 13*g12*g23**2/48 - 11*g13*g23**2/48)
    H_Trotter2 += sx1*sx3 * (13*g12**2*g13/48 - 11*g12**2*g23/48 - 1*g12*g13*g23/24 + 13*g13*g23**2/48 - 11*g12*g23**2/48)
    H_Trotter2 += sx2*sx3 * (13*g12**2*g23/48 - 11*g12**2*g13/48 - 1*g12*g13*g23/24 + 13*g13**2*g23/48 - 11*g12*g13**2/48)
    H_Trotter2 += sy1*sy2 * (1*g12*g13**2/12  - 1*g13**2*g23/24  - 1*g12*g13*g23/6  + 1*g12*g23**2/12  - 1*g13*g23**2/24 )
    H_Trotter2 += sy1*sy3 * (1*g12**2*g13/12  - 1*g12**2*g23/24  - 1*g12*g13*g23/6  + 1*g13*g23**2/12  - 1*g12*g23**2/24 )
    H_Trotter2 += sy2*sy3 * (1*g13**2*g23/12  - 1*g12**2*g13/24  - 1*g12*g13*g23/6  + 1*g12**2*g23/8   - 1*g12*g13**2/12 ) # different pattern
    H_Trotter2 += sz1*sz2 * (13*g13**2*g23/48 + 13*g13**2*g23/48 - 1*g12*g13*g23/24 - 11*g12*g13**2/48 - 11*g12*g23**2/48)
    H_Trotter2 += sz1*sz3 * (13*g12**2*g13/48 + 13*g12*g13**2/48 - 1*g12*g13*g23/24 - 11*g12**2*g23/48 - 11*g13**2*g23/48)

    # master equation solver
    res_ex = qt.mesolve(H_Heisenberg, psi0, tlist, [], [sz1, sz2, sz3])
    res_T1 = qt.mesolve(H_Trotter1, psi0, tlist, [], [sz1, sz2, sz3])
    res_T2 = qt.mesolve(H_Trotter2, psi0, tlist, [], [sz1, sz2, sz3])

    # expectation values
    evl_sz1_ex = res_ex.expect[0]
    evl_sz2_ex = res_ex.expect[1]
    evl_sz3_ex = res_ex.expect[2]
    evl_sz1_T1 = res_T1.expect[0]
    evl_sz2_T1 = res_T1.expect[1]
    evl_sz3_T1 = res_T1.expect[2]
    evl_sz1_T2 = res_T2.expect[0]
    evl_sz2_T2 = res_T2.expect[1]
    evl_sz3_T2 = res_T2.expect[2]

    inv_prob_1_ex = np.array([abs(evl_sz1_ex[0]-evl_sz1_ex[i])/2 for i in range(0, len(evl_sz1_ex))])
    inv_prob_2_ex = np.array([abs(evl_sz2_ex[0]-evl_sz2_ex[i])/2 for i in range(0, len(evl_sz2_ex))])
    inv_prob_3_ex = np.array([abs(evl_sz3_ex[0]-evl_sz3_ex[i])/2 for i in range(0, len(evl_sz3_ex))])
    inv_prob_1_T1 = np.array([abs(evl_sz1_T1[0]-evl_sz1_T1[i])/2 for i in range(0, len(evl_sz1_T1))])
    inv_prob_2_T1 = np.array([abs(evl_sz2_T1[0]-evl_sz2_T1[i])/2 for i in range(0, len(evl_sz2_T1))])
    inv_prob_3_T1 = np.array([abs(evl_sz3_T1[0]-evl_sz3_T1[i])/2 for i in range(0, len(evl_sz3_T1))])
    inv_prob_1_T2 = np.array([abs(evl_sz1_T2[0]-evl_sz1_T2[i])/2 for i in range(0, len(evl_sz1_T2))])
    inv_prob_2_T2 = np.array([abs(evl_sz2_T2[0]-evl_sz2_T2[i])/2 for i in range(0, len(evl_sz2_T2))])
    inv_prob_3_T2 = np.array([abs(evl_sz3_T2[0]-evl_sz3_T2[i])/2 for i in range(0, len(evl_sz3_T2))])

    # Plotting: Z expectation value (using exact model)
    fig, ax = plt.subplots(1, 1)
    fig.set_size_inches(16 / 2.54, 10 / 2.54)
    
    ax.plot(tlist, evl_sz1_ex, color="blue" , label=r"$\nu_1$")
    ax.plot(tlist, evl_sz2_ex, color="green", label=r"$\nu_2$")
    ax.plot(tlist, evl_sz3_ex, color="red"  , label=r"$\nu_3$")
    ax.set_title(r"Exact evolution of expectation value \langle Z_i\rangle")
    ax.set_xlabel("Time")
    ax.set_ylabel(r"Expectation value $\langle Z\rangle$")
    ax.set_xlim(np.min(tlist), np.max(tlist))
    legend = ax.legend(loc=3, frameon=True, borderaxespad=0.8, fontsize=8)
    legend.get_frame().set_facecolor('white')
    legend.get_frame().set_alpha(1.0)
    legend.get_frame().set_boxstyle("Square")

    plt.grid()
    plt.tight_layout()
    fig.savefig(f'..\\plots\\3q_flavor_inv\\Z_exp_vl.pdf', dpi=300)
    plt.close(fig)

    # Plotting: flavor inversion probability (using exact model)
    fig, ax = plt.subplots(1, 1)
    fig.set_size_inches(16 / 2.54, 10 / 2.54)
    
    ax.plot(tlist, inv_prob_1_ex, color="blue" , label=r"$\nu_1$")
    ax.plot(tlist, inv_prob_2_ex, color="green", label=r"$\nu_2$")
    ax.plot(tlist, inv_prob_3_ex, color="red"  , label=r"$\nu_3$")
    ax.set_title(r"Exact evolution of expectation value \langle Z_i\rangle")
    ax.set_xlabel("Time")
    ax.set_ylabel(r"Expectation value $\langle Z\rangle$")
    ax.set_xlim(np.min(tlist), np.max(tlist))
    legend = ax.legend(loc=3, frameon=True, borderaxespad=0.8, fontsize=8)
    legend.get_frame().set_facecolor('white')
    legend.get_frame().set_alpha(1.0)
    legend.get_frame().set_boxstyle("Square")

    plt.grid()
    plt.tight_layout()
    fig.savefig(f'..\\plots\\3q_flavor_inv\\track_flavor_inversion.pdf', dpi=300)
    plt.close(fig)
    




if __name__ == "__main__":
    main()