# ---------------------------- 3q_noisy_numop.py -----------------------------------
#-----------------------------------------------------------------------------------
# Compare evolution of the error in measuring an observable in a noisy scenario.
# I measure the occupancy number at different times and compare its expectation 
# value in various scenarios.


import numpy as np
import qutip as qt
import math
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors  
from matplotlib.offsetbox import AnchoredText
from iminuit import Minuit
from iminuit.cost import LeastSquares   
matplotlib.use('Agg')
import sys; sys.path.append("../classes")
from CommutatorTool import PauliCommutators


def main():
    ctl3 = PauliCommutators(N=3)
    T_tot = 20
    T_step = 0.2

    tlist = np.arange(0, T_tot, T_step)
    
    g12 = 1e-1
    g13 = 1e-1
    g23 = 2e-1

    sx1 = qt.tensor(qt.sigmax(), qt.qeye(2), qt.qeye(2))
    sy1 = qt.tensor(qt.sigmay(), qt.qeye(2), qt.qeye(2))
    sz1 = qt.tensor(qt.sigmaz(), qt.qeye(2), qt.qeye(2))
    sm1 = qt.tensor(qt.sigmam(), qt.qeye(2), qt.qeye(2))
    n1 = sm1 * sm1.dag()

    sx2 = qt.tensor(qt.qeye(2), qt.sigmax(), qt.qeye(2))
    sy2 = qt.tensor(qt.qeye(2), qt.sigmay(), qt.qeye(2))
    sz2 = qt.tensor(qt.qeye(2), qt.sigmaz(), qt.qeye(2))
    sm2 = qt.tensor(qt.qeye(2), qt.sigmam(), qt.qeye(2))
    n2 = sm2 * sm2.dag()

    sx3 = qt.tensor(qt.qeye(2), qt.qeye(2), qt.sigmax())
    sy3 = qt.tensor(qt.qeye(2), qt.qeye(2), qt.sigmay())
    sz3 = qt.tensor(qt.qeye(2), qt.qeye(2), qt.sigmaz())
    sm3 = qt.tensor(qt.qeye(2), qt.qeye(2), qt.sigmam())
    n3 = sm3 * sm3.dag()

    H_Heisenberg = g12*(sx1*sx2 + sy1*sy2 + sz1*sz2) + g13*(sx1*sx3 + sy1*sy3 + sz1*sz3) + g23*(sx2*sx3 + sy2*sy3 + sz2*sz3)

    H_Trotter1 = H_Heisenberg
    H_Trotter1 += 0.25 * (sx1*sz2*sy3 + sz1*sx2*sy3) * (g12*g13+g12*g23-3*g13*g23)
    H_Trotter1 += 0.25 * (sx1*sy2*sz3 + sz1*sy2*sz3) * (g12*g13+g13*g23-3*g12*g23)
    H_Trotter1 += 0.25 * (sy1*sx2*sz3 + sy1*sz2*sx3) * (g12*g23+g13*g23-3*g12*g13)

    psi0 = qt.tensor(qt.basis(2,1), qt.basis(2,1), (qt.basis(2,0) + qt.basis(2,1))/np.sqrt(2))

    H_Trotter2  = H_Heisenberg
    H_Trotter2 += sx1*sx2 * (13*g12*g13**2/48 - 11*g13**2*g23/48 - 1*g12*g13*g23/24 + 13*g12*g23**2/48 - 11*g13*g23**2/48)
    H_Trotter2 += sx1*sx3 * (13*g12**2*g13/48 - 11*g12**2*g23/48 - 1*g12*g13*g23/24 + 13*g13*g23**2/48 - 11*g12*g23**2/48)
    H_Trotter2 += sx2*sx3 * (13*g12**2*g23/48 - 11*g12**2*g13/48 - 1*g12*g13*g23/24 + 13*g13**2*g23/48 - 11*g12*g13**2/48)
    H_Trotter2 += sy1*sy2 * (1*g12*g13**2/12  - 1*g13**2*g23/24  - 1*g12*g13*g23/6  + 1*g12*g23**2/12  - 1*g13*g23**2/24 )
    H_Trotter2 += sy1*sy3 * (1*g12**2*g13/12  - 1*g12**2*g23/24  - 1*g12*g13*g23/6  + 1*g13*g23**2/12  - 1*g12*g23**2/24 )
    H_Trotter2 += sy2*sy3 * (1*g13**2*g23/12  - 1*g12**2*g13/24  - 1*g12*g13*g23/6  + 1*g12**2*g23/8   - 1*g12*g13**2/12 ) # different pattern
    H_Trotter2 += sz1*sz2 * (13*g13**2*g23/48 + 13*g13**2*g23/48 - 1*g12*g13*g23/24 - 11*g12*g13**2/48 - 11*g12*g23**2/48)
    H_Trotter2 += sz1*sz3 * (13*g12**2*g23/48 + 13*g12*g23**2/48 - 1*g12*g13*g23/24 - 11*g12**2*g13/48 - 11*g13*g23**2/48)
    H_Trotter2 += sz2*sz3 * (13*g12**2*g13/48 + 13*g12*g13**2/48 - 1*g12*g13*g23/24 - 11*g12**2*g23/48 - 11*g13**2*g23/48)

    # noise parameters
    # NOTE: T_relax < T_deph -> g_relax = 1 / T_relax > 1 / T_deph = g_deph
    g1_relax = np.sqrt(0.050)
    g1_deph  = np.sqrt(0.020)
    g2_relax = np.sqrt(0.030)
    g2_deph  = np.sqrt(0.020)
    g3_relax = np.sqrt(0.040)
    g3_deph  = np.sqrt(0.025)

    noise = [g1_relax*sm1.dag(), g2_relax*sm2.dag(), g3_relax*sm3.dag(), g1_deph*sz1, g2_deph*sz2, g3_deph*sz3]

    res_exact_HH = qt.mesolve(H_Heisenberg, psi0, tlist, [   ], [n1, n2, n3])
    res_noisy_HH = qt.mesolve(H_Heisenberg, psi0, tlist, noise, [n1, n2, n3])
    res_exact_T1 = qt.mesolve(H_Trotter1  , psi0, tlist, [   ], [n1, n2, n3])
    res_noisy_T1 = qt.mesolve(H_Trotter1  , psi0, tlist, noise, [n1, n2, n3])
    res_exact_T2 = qt.mesolve(H_Trotter2  , psi0, tlist, [   ], [n1, n2, n3])
    res_noisy_T2 = qt.mesolve(H_Trotter2  , psi0, tlist, noise, [n1, n2, n3])

    # plot: compare res_exact_HH with res_noisy_HH  
    fig, ax = plt.subplots(1, 1)
    fig.set_size_inches(16 / 2.54, 10 / 2.54)

    ax.plot(tlist, res_exact_HH.expect[0], color="#003CFF", label="Exact qubit 1")
    ax.plot(tlist, res_exact_HH.expect[1], color="#1EB600", label="Exact qubit 2")
    ax.plot(tlist, res_exact_HH.expect[2], color="#FF0000", label="Exact qubit 3")
    ax.plot(tlist, res_noisy_HH.expect[0], color="#00E0CE", linestyle="--", label="Noisy qubit 1")
    ax.plot(tlist, res_noisy_HH.expect[1], color="#00FF73", linestyle="--", label="Noisy qubit 2")
    ax.plot(tlist, res_noisy_HH.expect[2], color="#FF5050", linestyle="--", label="Noisy qubit 3")

    ax.set_title(r"Compare exact and noisy evolution of $\langle\sigma_-\sigma_+\rangle$ for $N=3$.")
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$\langle\sigma_-\sigma_+\rangle$")
    ax.set_xlim(np.min(tlist), np.max(tlist))

    legend = ax.legend(loc=1, frameon=True, borderaxespad=0.8, fontsize=8)
    legend.get_frame().set_facecolor('white')
    legend.get_frame().set_alpha(1.0)
    legend.get_frame().set_boxstyle("Square")

    plt.grid()
    plt.tight_layout()
    fig.savefig(f'..\\plots\\3q_noisy_numop\\compare_exact_and_noisy_Heisenberg.pdf', dpi=300)
    plt.close(fig)

    # plot: compare res_noisy_HH with res_noisy_T1
    fig, ax = plt.subplots(1, 1)
    fig.set_size_inches(16 / 2.54, 10 / 2.54)

    # ax.plot(tlist, res_noisy_HH.expect[0], color="#00E0CE", linestyle="--", label=r"Noisy HH $q_1$")
    # ax.plot(tlist, res_noisy_HH.expect[1], color="#00FF73", linestyle="--", label=r"Noisy HH $q_2$")
    # ax.plot(tlist, res_noisy_HH.expect[2], color="#FF5050", linestyle="--", label=r"Noisy HH $q_3$")
    # ax.plot(tlist, res_noisy_T1.expect[0], color="#007BE0", linestyle="-.", label=r"Noisy T1 $q_1$")
    # ax.plot(tlist, res_noisy_T1.expect[1], color="#88FF00", linestyle="-.", label=r"Noisy T1 $q_2$")
    # ax.plot(tlist, res_noisy_T1.expect[2], color="#FFA550", linestyle="-.", label=r"Noisy T1 $q_3$")

    # define Delta sigma_{A,B}^k = <s-1 s+1>_A - <s-1 s+1>_B
    ax.plot(tlist, res_noisy_HH.expect[0]-res_noisy_T1.expect[0], color="#3C00E0", linestyle="--", label=r"$\Delta\sigma_{H_\text{noisy},T1_\text{noisy}}^1$" )
    ax.plot(tlist, res_noisy_HH.expect[1]-res_noisy_T1.expect[1], color="#04C91E", linestyle="--", label=r"$\Delta\sigma_{H_\text{noisy},T1_\text{noisy}}^2$" )
    ax.plot(tlist, res_noisy_HH.expect[2]-res_noisy_T1.expect[2], color="#FF2600", linestyle="--", label=r"$\Delta\sigma_{H_\text{noisy},T1_\text{noisy}}^3$" )
    ax.plot(tlist, res_noisy_HH.expect[0]-res_noisy_T2.expect[0], color="#007BE0", linestyle="-.", label=r"$\Delta\sigma_{H_\text{noisy},T2_\text{noisy}}^1$" )
    ax.plot(tlist, res_noisy_HH.expect[1]-res_noisy_T2.expect[1], color="#88FF00", linestyle="-.", label=r"$\Delta\sigma_{H_\text{noisy},T2_\text{noisy}}^2$" )
    ax.plot(tlist, res_noisy_HH.expect[2]-res_noisy_T2.expect[2], color="#FFA550", linestyle="-.", label=r"$\Delta\sigma_{H_\text{noisy},T2_\text{noisy}}^3$" )
    ax.plot(tlist, res_noisy_T1.expect[0]-res_noisy_T2.expect[0], color="#00A0E0", linestyle=":" , label=r"$\Delta\sigma_{T1_\text{noisy},T2_\text{noisy}}^1$")
    ax.plot(tlist, res_noisy_T1.expect[1]-res_noisy_T2.expect[1], color="#00FF55", linestyle=":" , label=r"$\Delta\sigma_{T1_\text{noisy},T2_\text{noisy}}^2$")
    ax.plot(tlist, res_noisy_T1.expect[2]-res_noisy_T2.expect[2], color="#FF7950", linestyle=":" , label=r"$\Delta\sigma_{T1_\text{noisy},T2_\text{noisy}}^3$")

    ax.set_title(r"Compare noisy evolutions of $\langle\sigma_-\sigma_+\rangle$ for $N=3$.")
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$\Delta\sigma$")
    ax.set_xlim(np.min(tlist), np.max(tlist))

    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width, box.height])

    legend = ax.legend(loc=2, frameon=True, borderaxespad=0.8, fontsize=8, bbox_to_anchor=(1, 1))
    legend.get_frame().set_facecolor('white')
    legend.get_frame().set_alpha(1.0)
    legend.get_frame().set_boxstyle("Square")

    plt.grid()
    plt.tight_layout()
    fig.savefig(f'..\\plots\\3q_noisy_numop\\compare_noisy_H_T1_T2.pdf', dpi=300)
    plt.close(fig)



if __name__ == "__main__":
    main()
