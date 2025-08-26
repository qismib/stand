# --------------------------- 2q_ETH_Heisenberg.py ---------------------------
#-----------------------------------------------------------------------------
# Reproduce plots of https://arxiv.org/pdf/1502.06778 fig. 2(b), which
# displays the expectation values of the Pauli operators 
# sigma_{x,y,z}_{1,2} for the Heisenberg interaction as a function of the 
# quantum phase angle 2|J|τ along with the ideal evolution.
# In the paper, the dots indicate experimental data. I just perform a
# simulation using the exact model.

import numpy as np
import qutip as qt
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg')  



def main():
    # Parameters
    J = 1e-1 # interaction strength between qubit 1 and qubit 2
    max_phase_angle = 3*np.pi / 2
    tau_max = max_phase_angle / 2 / J

    tlist = np.linspace(0, tau_max, 500)

    # One qubit operators
    sx1 = qt.tensor(qt.sigmax(), qt.qeye(2))
    sy1 = qt.tensor(qt.sigmay(), qt.qeye(2))
    sz1 = qt.tensor(qt.sigmaz(), qt.qeye(2))

    sx2 = qt.tensor(qt.qeye(2), qt.sigmax())
    sy2 = qt.tensor(qt.qeye(2), qt.sigmay())
    sz2 = qt.tensor(qt.qeye(2), qt.sigmaz())

    # system Hamiltonian: isotropic Heisenberg interaction
    H_Heisenberg = J*(sx1*sx2 + sy1*sy2 + sz1*sz2)

    # initial state: |0> @ (|0> + |1>)/sqrt(2) (as stated in the paper)
    psi0 = qt.tensor(qt.basis(2,0), (qt.basis(2,0) + qt.basis(2,1)) / np.sqrt(2))

    res = qt.mesolve(H_Heisenberg, psi0, tlist, [], [sx1, sy1, sz1, sx2, sy2, sz2])
    ex1, ey1, ez1, ex2, ey2, ez2 = res.expect

    # Plotting
    fig, ax = plt.subplots(1, 3)
    fig.set_size_inches(16 / 2.54, 10 / 2.54)

    first = ax[0]
    secnd = ax[1]
    third = ax[2]

    phase_angle = 2*np.abs(J)*tlist # quantum phase angle (in [0, max_phase_angle])
    
    first.plot(phase_angle, ex1, 'r--', label=r"$\langle\sigma_x^1\rangle$")
    first.plot(phase_angle, ex2, 'b--', label=r"$\langle\sigma_x^2\rangle$")
    secnd.plot(phase_angle, ey1, 'r--', label=r"$\langle\sigma_y^1\rangle$") # not matching paper (overall -1 difference)
    secnd.plot(phase_angle, ey2, 'b--', label=r"$\langle\sigma_y^2\rangle$") # not matching paper (overall -1 difference)
    third.plot(phase_angle, ez1, 'r--', label=r"$\langle\sigma_z^1\rangle$")
    third.plot(phase_angle, ez2, 'b--', label=r"$\langle\sigma_z^2\rangle$")

    fig.suptitle("Expectation values of Pauli under Heisenberg.")

    first.set_xlabel(r'Phase angle $2|J|\tau$')
    first.set_xlim(np.min(phase_angle), np.max(phase_angle))
    first.set_ylim(-1, 1)
    
    secnd.set_xlabel(r'Phase angle $2|J|\tau$')
    secnd.set_xlim(np.min(phase_angle), np.max(phase_angle))
    secnd.set_ylim(-1, 1)

    third.set_xlabel(r'Phase angle $2|J|\tau$')
    third.set_xlim(np.min(phase_angle), np.max(phase_angle))
    third.set_ylim(-1, 1)

    first.set_ylabel("Expectation value.")

    x_ticks_labels = [r"0", r"$\frac{\pi}{2}$", r"$\pi$", r"$\frac{3\pi}{2}$"]
    y_ticks_labels = ["" for _ in range(1, 10)]

    first.set_xticks([0, np.pi/2, np.pi, max_phase_angle])
    secnd.set_xticks([0, np.pi/2, np.pi, max_phase_angle])
    third.set_xticks([0, np.pi/2, np.pi, max_phase_angle])
    first.set_xticklabels(x_ticks_labels, rotation='horizontal', fontsize=8)
    secnd.set_xticklabels(x_ticks_labels, rotation='horizontal', fontsize=8)
    third.set_xticklabels(x_ticks_labels, rotation='horizontal', fontsize=8)

    secnd.set_yticks(np.arange(-1, 1.25, 0.25))
    third.set_yticks(np.arange(-1, 1.25, 0.25))
    secnd.set_yticklabels(y_ticks_labels, rotation='horizontal', fontsize=2)
    third.set_yticklabels(y_ticks_labels, rotation='horizontal', fontsize=2)


    legend = first.legend(loc=3, frameon=True, borderaxespad=0.8, fontsize=8)
    legend.get_frame().set_facecolor('white')
    legend.get_frame().set_alpha(1.0)
    legend.get_frame().set_boxstyle("Square")

    legend = secnd.legend(loc=3, frameon=True, borderaxespad=0.8, fontsize=8)
    legend.get_frame().set_facecolor('white')
    legend.get_frame().set_alpha(1.0)
    legend.get_frame().set_boxstyle("Square")

    legend = third.legend(loc=3, frameon=True, borderaxespad=0.8, fontsize=8)
    legend.get_frame().set_facecolor('white')
    legend.get_frame().set_alpha(1.0)
    legend.get_frame().set_boxstyle("Square")

    first.grid(True, color="grey", linewidth="1.0", linestyle="--")
    secnd.grid(True, color="grey", linewidth="1.0", linestyle="--")
    third.grid(True, color="grey", linewidth="1.0", linestyle="--")

    plt.tight_layout()

    fig.savefig(f'..\\plots\\2q_ETH_Heisenberg\\expvl_Pauli.pdf', dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    main()