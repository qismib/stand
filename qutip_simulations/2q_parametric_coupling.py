# -------------------- 2q_num_parametric_coupling.py ------------------------
# ---------------------------------------------------------------------------
# Simulate a system of two qubit, one of which is tunable. I vary the 
# frequency of this qubit to 

from qutip import *
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
matplotlib.use('Agg')  


w1 = 4.5
w2 = 4

g = 0.1
DOmega = 1 * np.abs(w2 - w1) 

def w2_t(t, args=None):
    return w2 + DOmega * np.cos((w2 - w1)*t)


def main():
    t_list = np.linspace(-100, 100, 50000)

    # qubit 1 (ff)
    sz1 = tensor(sigmaz(), qeye(2))
    sx1 = tensor(sigmax(), qeye(2))
    sy1 = tensor(sigmay(), qeye(2))
    sm1 = tensor(sigmam(), qeye(2))
    n1 = sm1 * sm1.dag()

    # qubit 2 (tunable)
    sz2 = tensor(qeye(2), sigmaz())
    sx2 = tensor(qeye(2), sigmax())
    sy2 = tensor(qeye(2), sigmay())
    sm2 = tensor(qeye(2), sigmam())
    n2 = sm2 * sm2.dag()

    # Hamiltonian of the system
    Hconst = 0.5*w1*sz1 + g*(sx1 * sx2 + sy1 * sy2) # qubit 1 + interaction
    Ht = 0.5*sz2 # qubit 2
    H = [Hconst, [Ht, w2_t]]

    # initial state: 1 photon in qubit 1
    #                no photon in qubit 2
    psi0 = tensor(basis(2,1), basis(2,0)) # NOTE: IF CHANGE, ALSO CHANGE IN LEGEND OF fig1 PLOT

    res = mesolve(H, psi0, t_list, [], [n1, n2])

    # expectation values associated to number operators for the 2-dim subspace
    n1_res_list = res.expect[0]
    n2_res_list = res.expect[1]

    # --------------------------------------------------------------------------------------------------------------------
    # plot exp_vl vs. t   

    # Create figure and subplots using gridspec
    fig1, axis = plt.subplots(2, 1, gridspec_kw={"height_ratios": [3, 1]})
    fig1.set_size_inches(16 / 2.54, 10 / 2.54)
    first = axis[0]
    secnd = axis[1]

    # Plot area
    first.plot(t_list, n1_res_list, color="blue", label=r"$\langle \sigma_-^1\sigma_+^1\rangle(t)$")
    first.plot(t_list, n2_res_list, color="red" , label=r"$\langle \sigma_-^2\sigma_+^2\rangle(t)$")
    first.set_title("Time evolution of number operator")
    first.set_xlabel("Time")
    first.set_ylabel(r'$\langle \sigma_-\sigma_+\rangle$')
    first.set_xlim(np.min(t_list), np.max(t_list))
    first.grid(linestyle='--', linewidth=0.5)

    legend = first.legend(loc=3, frameon=True, borderaxespad=0.8, fontsize=8)
    legend.get_frame().set_facecolor('white')
    legend.get_frame().set_edgecolor('black')
    legend.get_frame().set_alpha(1.0)
    legend.get_frame().set_boxstyle("Square")

    # Table area
    secnd.axis('off')  # Hide axes

    # Custom text content for the table
    text = [
        [f"$\\omega_1 = {w1}$, $\\omega_2 \\equiv \\omega_2 (t = 0) = {w2}$, $\\Delta(t=0) = {np.abs(w2 - w1)}$, $g = {g}$, $\\Delta\\omega = {DOmega}$."],
        ["$\\omega_2 (t) = \\omega_2 + \\Delta\\omega \\dot \\cos((\\omega_2 - \\omega_1)t)$"],
        ["Initial state: $|\\psi_0\\rangle = |1\\rangle \\otimes |0\\rangle$"]
    ]

    # Create table with 1 column
    table = secnd.table(cellText=text, colLabels=["Parameters of the simulation"], loc="center", cellLoc="left")

    # Optional styling
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1.2, 1.2)  # scale width and height

    # Adjust layout
    plt.tight_layout()
    fig1.savefig(f'..\\plots\\2q_parametric_coupling\\n1n2_vs_time_w1_{w1}_w2_{w2}_g_{g}_Delta_{np.abs(w2 - w1)}_DOmega_{DOmega}_v2.pdf', dpi=300) # Delta is |w2 - w1|(t=0)


    

if __name__ == "__main__":
    main()