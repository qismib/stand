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

    tlist = np.linspace(0, 20, 1000)

    sx1 = qt.tensor(qt.sigmax(), qt.qeye(2), qt.qeye(2))
    sy1 = qt.tensor(qt.sigmay(), qt.qeye(2), qt.qeye(2))
    sz1 = qt.tensor(qt.sigmaz(), qt.qeye(2), qt.qeye(2))

    sx2 = qt.tensor(qt.qeye(2), qt.sigmax(), qt.qeye(2))
    sy2 = qt.tensor(qt.qeye(2), qt.sigmay(), qt.qeye(2))
    sz2 = qt.tensor(qt.qeye(2), qt.sigmaz(), qt.qeye(2))

    sx3 = qt.tensor(qt.qeye(2), qt.qeye(2), qt.sigmax())
    sy3 = qt.tensor(qt.qeye(2), qt.qeye(2), qt.sigmay())
    sz3 = qt.tensor(qt.qeye(2), qt.qeye(2), qt.sigmaz())

    E_Trotter1  = 0.25 * (sx1*sz2*sy3 + sz1*sx2*sy3) * (g12*g13+g12*g23-3*g13*g23)
    E_Trotter1 += 0.25 * (sx1*sy2*sz3 + sz1*sy2*sz3) * (g12*g13+g13*g23-3*g12*g23)
    E_Trotter1 += 0.25 * (sy1*sx2*sz3 + sy1*sz2*sx3) * (g12*g23+g13*g23-3*g12*g13)

    E_Trotter2  = sx1*sx2 * (13*g12*g13**2/48 - 11*g13**2*g23/48 - 1*g12*g13*g23/24 + 13*g12*g23**2/48 - 11*g13*g23**2/48)
    E_Trotter2 += sx1*sx3 * (13*g12**2*g13/48 - 11*g12**2*g23/48 - 1*g12*g13*g23/24 + 13*g13*g23**2/48 - 11*g12*g23**2/48)
    E_Trotter2 += sx2*sx3 * (13*g12**2*g23/48 - 11*g12**2*g13/48 - 1*g12*g13*g23/24 + 13*g13**2*g23/48 - 11*g12*g13**2/48)
    E_Trotter2 += sy1*sy2 * (1*g12*g13**2/12  - 1*g13**2*g23/24  - 1*g12*g13*g23/6  + 1*g12*g23**2/12  - 1*g13*g23**2/24 )
    E_Trotter2 += sy1*sy3 * (1*g12**2*g13/12  - 1*g12**2*g23/24  - 1*g12*g13*g23/6  + 1*g13*g23**2/12  - 1*g12*g23**2/24 )
    E_Trotter2 += sy2*sy3 * (1*g13**2*g23/12  - 1*g12**2*g13/24  - 1*g12*g13*g23/6  + 1*g12**2*g23/8   - 1*g12*g13**2/12 ) # different pattern
    E_Trotter2 += sz1*sz2 * (13*g13**2*g23/48 + 13*g13**2*g23/48 - 1*g12*g13*g23/24 - 11*g12*g13**2/48 - 11*g12*g23**2/48)
    E_Trotter2 += sz1*sz3 * (13*g12**2*g13/48 + 13*g12*g13**2/48 - 1*g12*g13*g23/24 - 11*g12**2*g23/48 - 11*g13**2*g23/48)

    H_Heisenberg = g12*(sx1*sx2 + sy1*sy2 + sz1*sz2) + g13*(sx1*sx3 + sy1*sy3 + sz1*sz3) + g23*(sx2*sx3 + sy2*sy3 + sz2*sz3)

    normA_Trotter1 = np.array([np.max(np.linalg.svd(((-1j*(H_Heisenberg + E_Trotter1*t)*t).expm() - (-1j*H_Heisenberg*t).expm()).full(), compute_uv=False)) for t in tlist])
    normA_Trotter2 = np.array([np.max(np.linalg.svd(((-1j*(H_Heisenberg + E_Trotter2*t*t)*t).expm() - (-1j*H_Heisenberg*t).expm()).full(), compute_uv=False)) for t in tlist])

    fig, ax = plt.subplots(1, 1)
    fig.set_size_inches(16 / 2.54, 10 / 2.54)

    ax.plot(tlist, normA_Trotter1, color="#093DE9", label="First order Trotter")
    ax.plot(tlist, normA_Trotter2, color="#031E77", label="Second order Trotter")
    ax.set_title(r"Compare evolution of additive Trotter error.")
    ax.set_xlabel("Time")
    ax.set_ylabel(r"$\|e^{-iH_{\text{Trotter}}t} - e^{-iH_\text{Heisenberg}t}\|$")
    ax.set_xlim(np.min(tlist), np.max(tlist))
    legend = ax.legend(loc=3, frameon=True, borderaxespad=0.8, fontsize=8)
    legend.get_frame().set_facecolor('white')
    legend.get_frame().set_alpha(1.0)
    legend.get_frame().set_boxstyle("Square")

    plt.grid()
    plt.tight_layout()
    fig.savefig(f'..\\plots\\3q_Trotter_compare\\additive_err_evo.pdf', dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    main()