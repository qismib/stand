# ------------------------- 3q_evo_Trotter1.py -------------------------
#----------------------------------------------------------------------
# Compare the propagators U = e^{-iHt} with for H_Heisenberg and
# the first order Trotter Hamiltonian, which contains spurios terms
# in addition to the correct Heisenberg contriutions.

import numpy as np
import qutip as qt
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg')  




def main():
    g12 = 1e-1
    g13 = 2e-1
    g23 = 1e-1

    tlist = np.linspace(0, 10, 1000)

    sx1 = qt.tensor(qt.sigmax(), qt.qeye(2), qt.qeye(2))
    sy1 = qt.tensor(qt.sigmay(), qt.qeye(2), qt.qeye(2))
    sz1 = qt.tensor(qt.sigmaz(), qt.qeye(2), qt.qeye(2))

    sx2 = qt.tensor(qt.qeye(2), qt.sigmax(), qt.qeye(2))
    sy2 = qt.tensor(qt.qeye(2), qt.sigmay(), qt.qeye(2))
    sz2 = qt.tensor(qt.qeye(2), qt.sigmaz(), qt.qeye(2))

    sx3 = qt.tensor(qt.qeye(2), qt.qeye(2), qt.sigmax())
    sy3 = qt.tensor(qt.qeye(2), qt.qeye(2), qt.sigmay())
    sz3 = qt.tensor(qt.qeye(2), qt.qeye(2), qt.sigmaz())

    psi0 = qt.tensor(qt.basis(2,0), qt.basis(2,0), qt.basis(2,0)) # |000>

    E_Trotter2  = sx1*sx2 * (13*g12*g13**2/48 - 11*g13**2*g23/48 - 1*g12*g13*g23/24 + 13*g12*g23**2/48 - 11*g13*g23**2/48)
    E_Trotter2 += sx1*sx3 * (13*g12**2*g13/48 - 11*g12**2*g23/48 - 1*g12*g13*g23/24 + 13*g13*g23**2/48 - 11*g12*g23**2/48)
    E_Trotter2 += sx2*sx3 * (13*g12**2*g23/48 - 11*g12**2*g13/48 - 1*g12*g13*g23/24 + 13*g13**2*g23/48 - 11*g12*g13**2/48)
    E_Trotter2 += sy1*sy2 * (1*g12*g13**2/12  - 1*g13**2*g23/24  - 1*g12*g13*g23/6  + 1*g12*g23**2/12  - 1*g13*g23**2/24 )
    E_Trotter2 += sy1*sy3 * (1*g12**2*g13/12  - 1*g12**2*g23/24  - 1*g12*g13*g23/6  + 1*g13*g23**2/12  - 1*g12*g23**2/24 )
    E_Trotter2 += sy2*sy3 * (1*g13**2*g23/12  - 1*g12**2*g13/24  - 1*g12*g13*g23/6  + 1*g12**2*g23/8   - 1*g12*g13**2/12 ) # different pattern
    E_Trotter2 += sz1*sz2 * (13*g13**2*g23/48 + 13*g13**2*g23/48 - 1*g12*g13*g23/24 - 11*g12*g13**2/48 - 11*g12*g23**2/48)
    E_Trotter2 += sz1*sz3 * (13*g12**2*g13/48 + 13*g12*g13**2/48 - 1*g12*g13*g23/24 - 11*g12**2*g23/48 - 11*g13**2*g23/48)


    H_Heisenberg = g12*(sx1*sx2 + sy1*sy2 + sz1*sz2) + g13*(sx1*sx3 + sy1*sy3 + sz1*sz3) + g23*(sx2*sx3 + sy2*sy3 + sz2*sz3)

    # #1: compute expectation value of the two propagators over some fixed state
    U_Heisenberg = np.array([qt.expect((-1j*H_Heisenberg*t).expm(), psi0).real for t in tlist])
    U_Trotter2   = np.array([qt.expect((-1j*(H_Heisenberg + E_Trotter2*t)*t).expm(), psi0).real for t in tlist])

    fig, ax = plt.subplots(1, 1)
    fig.set_size_inches(16 / 2.54, 10 / 2.54)

    ax.plot(tlist, U_Heisenberg, color="blue", label=r"$\exp{-iH_{Heisenberg}t}$")
    ax.plot(tlist, U_Trotter2, color="green", label=r"$\exp{-iH_{\text{Trotter }2}t}$")
    ax.set_title(r"Evolution of exact and second order Trotter propagator.")
    ax.set_xlabel("Time")
    ax.set_ylabel(r"$\langle e^{-iHt}\rangle$")
    ax.set_xlim(np.min(tlist), np.max(tlist))
    # legend = ax.legend(loc=3, frameon=True, borderaxespad=0.8, fontsize=8)
    # legend.get_frame().set_facecolor('white')
    # legend.get_frame().set_alpha(1.0)
    # legend.get_frame().set_boxstyle("Square")

    plt.grid()
    plt.tight_layout()
    fig.savefig(f'..\\plots\\3q_evo_Trotter2\\comp_exp_val.pdf', dpi=300)
    plt.close(fig)

    # #2: compute the spectral norm difference ||A(t)||=||e^{-iH_Trotter t}-e^{-iH_Heisenberg t}||

    normA_Trotter1 = [np.sqrt(abs(np.max((((-1j*(H_Heisenberg + E_Trotter2*t)*t).expm() - (-1j*H_Heisenberg*t).expm()).eigenstates()[0])))) for t in tlist]

    fig, ax = plt.subplots(1, 1)
    fig.set_size_inches(16 / 2.54, 10 / 2.54)

    ax.plot(tlist, normA_Trotter1, color="blue")
    ax.set_title(r"Evolution of first order additive Trotter error.")
    ax.set_xlabel("Time")
    ax.set_ylabel(r"$\|e^{-iH_{\text{Trotter }2}t} - e^{-iH_\text{Heisenberg}t}\|$")
    ax.set_xlim(np.min(tlist), np.max(tlist))
    # legend = ax.legend(loc=3, frameon=True, borderaxespad=0.8, fontsize=8)
    # legend.get_frame().set_facecolor('white')
    # legend.get_frame().set_alpha(1.0)
    # legend.get_frame().set_boxstyle("Square")

    plt.grid()
    plt.tight_layout()
    fig.savefig(f'..\\plots\\3q_evo_Trotter2\\additive_err_evo.pdf', dpi=300)
    plt.close(fig)

    # #3: fidelity approach F(t) = ||<psi0|U_Trotter^\dag * U_Heisenberg|psi0>||

    fidelity_val = np.zeros(len(tlist))

    for i, t in enumerate(tlist):
        U_Heis = (-1j*H_Heisenberg*t).expm()
        U_Tr1  = (-1j*(H_Heisenberg + E_Trotter2*t)*t).expm()
        fidelity_val[i] = abs((U_Tr1*psi0).dag() * (U_Heis*psi0))

    fig, ax = plt.subplots(1, 1)
    fig.set_size_inches(16 / 2.54, 10 / 2.54)

    ax.plot(tlist, fidelity_val, color="blue")
    ax.set_title(r"Evolution of fidelity.")
    ax.set_xlabel("Time")
    ax.set_ylabel(r"$\mathcal{F}(t) = \|\langle\psi_0\vert U_{\text{Trotter }2}^\dag - U_\text{Heisenberg}\vert\psi_0\rangle\|$")
    ax.set_xlim(np.min(tlist), np.max(tlist))
    # legend = ax.legend(loc=3, frameon=True, borderaxespad=0.8, fontsize=8)
    # legend.get_frame().set_facecolor('white')
    # legend.get_frame().set_alpha(1.0)
    # legend.get_frame().set_boxstyle("Square")

    plt.grid()
    plt.tight_layout()
    fig.savefig(f'..\\plots\\3q_evo_Trotter2\\fidelity_evo.pdf', dpi=300)
    plt.close(fig)



if __name__ == "__main__":
    main()