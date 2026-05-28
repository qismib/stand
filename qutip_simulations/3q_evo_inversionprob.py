# ------------------------- 3q_evo_inversionprob.py -------------------------
#----------------------------------------------------------------------------
# Simulate a system of 3 neutrinos and track the inversion probability as
# a function of time. I compare the exact Heisenberg Hamiltonian with the 
# first and second order Trotter Hamilonians.

import numpy as np
import qutip as qt
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg')  






def main():

    g12 = 1e-1
    g13 = 2e-1
    g23 = 1e-1

    tlist = np.linspace(0, 100, 10000)

    sx1 = qt.tensor(qt.sigmax(), qt.qeye(2), qt.qeye(2))
    sy1 = qt.tensor(qt.sigmay(), qt.qeye(2), qt.qeye(2))
    sz1 = qt.tensor(qt.sigmaz(), qt.qeye(2), qt.qeye(2))

    sx2 = qt.tensor(qt.qeye(2), qt.sigmax(), qt.qeye(2))
    sy2 = qt.tensor(qt.qeye(2), qt.sigmay(), qt.qeye(2))
    sz2 = qt.tensor(qt.qeye(2), qt.sigmaz(), qt.qeye(2))

    sx3 = qt.tensor(qt.qeye(2), qt.qeye(2), qt.sigmax())
    sy3 = qt.tensor(qt.qeye(2), qt.qeye(2), qt.sigmay())
    sz3 = qt.tensor(qt.qeye(2), qt.qeye(2), qt.sigmaz())

    
    H_Heisenberg = g12*(sx1*sx2 + sy1*sy2 + sz1*sz2) + g13*(sx1*sx3 + sy1*sy3 + sz1*sz3) + g23*(sx2*sx3 + sy2*sy3 + sz2*sz3)

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

    psi0 = qt.tensor(qt.basis(2,0), qt.basis(2,1), qt.basis(2,0)) # |010>

    # I choose to track only the evolution of neutrino 1
    res_ex = qt.mesolve(H_Heisenberg, psi0, tlist, [], [sz1])
    res_T1 = qt.mesolve(H_Trotter1, psi0, tlist, [], [sz1])
    res_T2 = qt.mesolve(H_Trotter2, psi0, tlist, [], [sz1])

    evl_sz_ex = res_ex.expect[0]
    evl_sz_T1 = res_T1.expect[0]
    evl_sz_T2 = res_T2.expect[0]

    # Flavor inversion probability
    # equation (13) in https://arxiv.org/pdf/2207.03189 defines the inversion probability as P_inv(t) = |<Z>(0) - <Z>(t)|/2, where <Z>(t) = <psi(t)|Z|psi(t)>
    inv_prob_ex = np.array([abs(evl_sz_ex[0]-evl_sz_ex[i])/2 for i in range(0, len(evl_sz_ex))])
    inv_prob_T1 = np.array([abs(evl_sz_T1[0]-evl_sz_T1[i])/2 for i in range(0, len(evl_sz_T1))])
    inv_prob_T2 = np.array([abs(evl_sz_T2[0]-evl_sz_T2[i])/2 for i in range(0, len(evl_sz_T2))])


    # Plots
    fig, ax = plt.subplots(1, 1)
    fig.set_size_inches(16 / 2.54, 10 / 2.54)

    ax.plot(tlist, inv_prob_ex, color="blue", label="Exact Heisenberg")
    ax.plot(tlist, inv_prob_T1, color="green", label="First order Trotter")
    ax.plot(tlist, inv_prob_T2, color="red", label="Second order Trotter")
    ax.set_title(r"Flavor inversion probability.")
    ax.set_xlabel("Time")
    ax.set_ylabel(r"$P_\text{inv}(t)=\frac{|Z(0)-Z(t)|}{2}$")
    ax.set_xlim(np.min(tlist), np.max(tlist))
    legend = ax.legend(loc=3, frameon=True, borderaxespad=0.8, fontsize=8)
    legend.get_frame().set_facecolor('white')
    legend.get_frame().set_alpha(1.0)
    legend.get_frame().set_boxstyle("Square")

    plt.grid()
    plt.tight_layout()
    fig.savefig(f'..\\plots\\3q_evo_inversionprob\\track_flavor_inversion.pdf', dpi=300)
    plt.close(fig)
    




if __name__ == "__main__":
    main()