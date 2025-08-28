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

    E_Trotter1  = 0.25 * (sx1*sz2*sy3 + sz1*sx2*sy3) * (g12*g13+g12*g23-3*g13*g23)
    E_Trotter1 += 0.25 * (sx1*sy2*sz3 + sz1*sy2*sz3) * (g12*g13+g13*g23-3*g12*g23)
    E_Trotter1 += 0.25 * (sy1*sx2*sz3 + sy1*sz2*sx3) * (g12*g23+g13*g23-3*g12*g13)

    H_Heisenberg = g12*(sx1*sx2 + sy1*sy2 + sz1*sz2) + g13*(sx1*sx3 + sy1*sy3 + sz1*sz3) + g23*(sx2*sx3 + sy2*sy3 + sz2*sz3)

    # #1: compute expectation value of the two propagators over some fixed state
    U_Heisenberg = np.array([qt.expect((-1j*H_Heisenberg*t).expm(), psi0).real for t in tlist])
    U_Trotter1   = np.array([qt.expect((-1j*(H_Heisenberg + E_Trotter1*t)*t).expm(), psi0).real for t in tlist])

    fig, ax = plt.subplots(1, 1)
    fig.set_size_inches(16 / 2.54, 10 / 2.54)

    ax.plot(tlist, U_Heisenberg, color="blue", label=r"$\exp{-iH_{Heisenberg}t}$")
    ax.plot(tlist, U_Trotter1, color="green", label=r"$\exp{-iH_{\text{Trotter }1}t}$")
    ax.set_title(r"Evolution of exact and first order Trotter propagator.")
    ax.set_xlabel("Time")
    ax.set_ylabel(r"$e^{-iHt}$")
    ax.set_xlim(np.min(tlist), np.max(tlist))
    # legend = ax.legend(loc=3, frameon=True, borderaxespad=0.8, fontsize=8)
    # legend.get_frame().set_facecolor('white')
    # legend.get_frame().set_alpha(1.0)
    # legend.get_frame().set_boxstyle("Square")

    plt.grid()
    plt.tight_layout()
    fig.savefig(f'..\\plots\\3q_evo_Trotter1\\comp_exp_val.pdf', dpi=300)
    plt.close(fig)

    # #2: plot the spectral norm difference ||A||(t)=||e^{-iH_Trotter t}-e^{-iH_Heisenberg t}|| for t in tlist

    normA_Trotter1 = np.array([np.max(np.linalg.svd(((-1j*(H_Heisenberg + E_Trotter1*t)*t).expm() - (-1j*H_Heisenberg*t).expm()).full(), compute_uv=False)) for t in tlist])

    fig, ax = plt.subplots(1, 1)
    fig.set_size_inches(16 / 2.54, 10 / 2.54)

    ax.plot(tlist, normA_Trotter1, color="blue")
    ax.set_title(r"Evolution of first order additive Trotter error.")
    ax.set_xlabel("Time")
    ax.set_ylabel(r"$\|e^{-iH_{\text{Trotter }1}t} - e^{-iH_\text{Heisenberg}t}\|$")
    ax.set_xlim(np.min(tlist), np.max(tlist))
    # legend = ax.legend(loc=3, frameon=True, borderaxespad=0.8, fontsize=8)
    # legend.get_frame().set_facecolor('white')
    # legend.get_frame().set_alpha(1.0)
    # legend.get_frame().set_boxstyle("Square")

    plt.grid()
    plt.tight_layout()
    fig.savefig(f'..\\plots\\3q_evo_Trotter1\\additive_err_evo.pdf', dpi=300)
    plt.close(fig)

    # #3: fidelity approach F(t) = ||<psi0|U_Trotter^\dag * U_Heisenberg|psi0>||

    fidelity_val = np.zeros(len(tlist))

    for i, t in enumerate(tlist):
        U_Heis = (-1j*H_Heisenberg*t).expm()
        U_Tr1  = (-1j*(H_Heisenberg + E_Trotter1*t)*t).expm()
        fidelity_val[i] = abs((U_Tr1*psi0).dag() * (U_Heis*psi0))

    fig, ax = plt.subplots(1, 1)
    fig.set_size_inches(16 / 2.54, 10 / 2.54)

    ax.plot(tlist, fidelity_val, color="blue")
    ax.set_title(r"Evolution of fidelity.")
    ax.set_xlabel("Time")
    ax.set_ylabel(r"$\mathcal{F}(t) = \|\langle\psi_0\vert U_{\text{Trotter }1}^\dag - U_\text{Heisenberg}\vert\psi_0\rangle\|$")
    ax.set_xlim(np.min(tlist), np.max(tlist))
    # legend = ax.legend(loc=3, frameon=True, borderaxespad=0.8, fontsize=8)
    # legend.get_frame().set_facecolor('white')
    # legend.get_frame().set_alpha(1.0)
    # legend.get_frame().set_boxstyle("Square")

    plt.grid()
    plt.tight_layout()
    fig.savefig(f'..\\plots\\3q_evo_Trotter1\\fidelity_evo.pdf', dpi=300)
    plt.close(fig)



if __name__ == "__main__":
    main()