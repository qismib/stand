# ------------------------- 2q_num_op_comparison.py -------------------------
# ---------------------------------------------------------------------------
# Using the XY evolution as a toy model, I compare the evolution of the
# number operator for a two-dimensional system in different scenarios.


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
    sm_man1 = tensor(0.5*(sigmax() - 1j*sigmay()), qeye(2)) # "manually" defined sigma-
    sm1 = tensor(sigmam(), qeye(2)) # qutip built in sigma-
    n1 = sm1 * sm1.dag()
    n1_man = sm_man1 * sm_man1.dag()

    # qubit 2 (tunable)
    sz2 = tensor(qeye(2), sigmaz())
    sm2 = tensor(qeye(2), sigmam())
    sx2 = tensor(qeye(2), sigmax())
    sy2 = tensor(qeye(2), sigmay())
    sm_man2 = tensor(qeye(2), 0.5*(sigmax() - 1j*sigmay())) # "manually" defined sigma-
    n2 = sm2 * sm2.dag() # qutip built in sigma-
    n2_man = sm_man2 * sm_man2.dag()

    # Hamiltonian of the system
    Hs = 0.5*w1*sz1 + g*(sm1.dag()*sm2 + sm1*sm2.dag()) # qubit 1 + s+s- interaction 
    Hg = 0.5*w1*sz1 + 0.5*g*(sx1*sx2 + sy1*sy2) # qubit 1 + XY interaction
    Ht = 0.5*sz2 # qubit 2
    H1 = [Hs, [Ht, w2_t]]
    H2 = [Hg, [Ht, w2_t]]

    # initial state: 1 photon in qubit 1
    #                0 photon in qubit 2
    psi0 = tensor(basis(2,1), basis(2,0)) # NOTE: IF CHANGE, ALSO CHANGE IN LEGEND OF fig1 PLOT

    # evolution according to sigma+- interaction hamiltonian
    res1 = mesolve(H1, psi0, t_list, [], [n1, n2, n1_man, n2_man])

    # evolution according to XY interaction hamiltonian
    res2 = mesolve(H2, psi0, t_list, [], [n1, n2, n1_man, n2_man])

    n1_res1_list = res1.expect[0]
    n2_res1_list = res1.expect[1]
    n1_man_res1_list = res1.expect[2]
    n2_man_res1_list = res1.expect[3]

    n1_res2_list = res2.expect[0]
    n2_res2_list = res2.expect[1]
    n1_man_res2_list = res2.expect[2]
    n2_man_res2_list = res2.expect[3]

    # --------------------------------------------------------------------------------------------------------------------
    # plot exp_vl qubit 1 vs. t   

    fig1, ax = plt.subplots(2, 2)
    fig1.set_size_inches(16 / 2.54, 10 / 2.54)

    oneone = ax[0][0]
    onetwo = ax[0][1]
    twoone = ax[1][0]
    twotwo = ax[1][1]
    
    oneone.plot(t_list, n1_res1_list, color="blue")
    onetwo.plot(t_list, n1_man_res1_list, color="cyan")
    twoone.plot(t_list, n1_res2_list, color="red")
    twotwo.plot(t_list, n1_man_res2_list, color="orange")

    oneone.set_title(r"$\langle n_1(t)\rangle$ with $H_{XY}$")
    onetwo.set_title(r"$\langle n_1^\text{man}(t)\rangle$ with $H_{XY}$")
    twoone.set_title(r"$\langle n_1(t)\rangle$ with $H_{\sigma_\pm}$")
    twotwo.set_title(r"$\langle n_1^\text{man}(t)\rangle$ with $H_{\sigma_\pm}$")

    oneone.set_xlabel("Time")
    onetwo.set_xlabel("Time")
    twoone.set_xlabel("Time")
    twotwo.set_xlabel("Time")
    oneone.set_ylabel(r"$\langle n_1(t)\rangle$")
    onetwo.set_ylabel(r"$\langle n_1^\text{man}(t)\rangle$")
    twoone.set_ylabel(r"$\langle n_1(t)\rangle$")
    twotwo.set_ylabel(r"$\langle n_1^\text{man}(t)\rangle$")

    oneone.set_xlim(np.min(t_list), np.max(t_list))
    onetwo.set_xlim(np.min(t_list), np.max(t_list))
    twoone.set_xlim(np.min(t_list), np.max(t_list))
    twotwo.set_xlim(np.min(t_list), np.max(t_list))

    oneone.grid(linestyle='--', linewidth=0.5)
    onetwo.grid(linestyle='--', linewidth=0.5)
    twoone.grid(linestyle='--', linewidth=0.5)
    twotwo.grid(linestyle='--', linewidth=0.5)

    # legoneone = oneone.legend(loc=3, frameon=True, borderaxespad=0.8, fontsize=8)
    # legonetwo = onetwo.legend(loc=3, frameon=True, borderaxespad=0.8, fontsize=8)
    # legtwoone = twoone.legend(loc=3, frameon=True, borderaxespad=0.8, fontsize=8)
    # legtwotwo = twotwo.legend(loc=3, frameon=True, borderaxespad=0.8, fontsize=8)
    # legoneone.get_frame().set_facecolor('white')
    # legonetwo.get_frame().set_facecolor('white')
    # legtwoone.get_frame().set_facecolor('white')
    # legtwotwo.get_frame().set_facecolor('white')
    # legoneone.get_frame().set_edgecolor('black')
    # legonetwo.get_frame().set_edgecolor('black')
    # legtwoone.get_frame().set_edgecolor('black')
    # legtwotwo.get_frame().set_edgecolor('black')
    # legoneone.get_frame().set_alpha(1.0)
    # legonetwo.get_frame().set_alpha(1.0)
    # legtwoone.get_frame().set_alpha(1.0)
    # legtwotwo.get_frame().set_alpha(1.0)
    # legoneone.get_frame().set_boxstyle("Square")
    # legonetwo.get_frame().set_boxstyle("Square")
    # legtwoone.get_frame().set_boxstyle("Square")
    # legtwotwo.get_frame().set_boxstyle("Square")

    plt.tight_layout()
    fig1.savefig(f'..\\plots\\2q_num_op_comparison\\comparison_number_ops.pdf', dpi=300)

if __name__ == "__main__":
    main()