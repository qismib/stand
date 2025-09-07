# ------------------------- 3q_Trotter1_compare_bounds.py -------------------------
#----------------------------------------------------------------------------------
# Compare bounds for the additive first order Trotter error obtained in section 
# 4.1.1 with a qutip simulation of the system evolving under the first order
# propagator.

import numpy as np
import qutip as qt
import sys; sys.path.append("../classes")
from iminuit import Minuit
from iminuit.cost import LeastSquares
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
matplotlib.use('Agg')  




def alpha_1(t, g12, g13, g23):
    coeff  = 2 * (abs(g12*g13-g12*g23+g13*g23) + abs(g12*g13+g12*g23-g13*g23 ) + abs(-g12*g13+g12*g23+g13*g23))
    coeff += 4 * (abs(g12*g13-g12*g23-g13*g23) + abs(-g12*g13-g12*g23+g13*g23) + abs(-g12*g13+g12*g23-g13*g23))
    return coeff * t**2

def alpha_2(t, g12, g13, g23):
    return 18 * (g12*g13 + g12*g23 + g13*g23) * t**2

def alpha_3(t, g12, g13, g23):
    eta = 3.464 # TODO integrate exact result obtained in other script
    return 3*eta * (g12*g13 + g12*g23 + g13*g23) * t**2

def power_law(x, beta, gamma): # y = beta * x^gamma
    return beta * x**gamma

# ------------------- MAIN -------------------

def main():
    T = 0.21
    r = 50

    g12 = 1e-1
    g13 = 1e-1
    g23 = 2e-1

    tlist = np.linspace(0.01, T, r)

    # simulation of a system of N=3 qubits evolving according to the Trotterized first order propagator
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

    E_Trotter1  = 0.25 * (sx1*sz2*sy3 + sz1*sx2*sy3) * (g12*g13+g12*g23-3*g13*g23)
    E_Trotter1 += 0.25 * (sx1*sy2*sz3 + sz1*sy2*sz3) * (g12*g13+g13*g23-3*g12*g23)
    E_Trotter1 += 0.25 * (sy1*sx2*sz3 + sy1*sz2*sx3) * (g12*g23+g13*g23-3*g12*g13)

    # compute norm distances between exact and Trotterized propagators

    normA_Trotter1 = np.array([np.max(np.linalg.svd(((-1j*(H_Heisenberg + E_Trotter1*t)*t).expm() - (-1j*H_Heisenberg*t).expm()).full(), compute_uv=False)) for t in tlist])

    fig, ax = plt.subplots(1, 1)
    fig.set_size_inches(16 / 2.54, 10 / 2.54)

    # fit data
    alpha1_list = alpha_1(tlist, g12, g13, g23)
    alpha2_list = alpha_2(tlist, g12, g13, g23)
    alpha3_list = alpha_3(tlist, g12, g13, g23)

    ls_alpha1 = LeastSquares(tlist, alpha1_list, 0*tlist + 0.01, power_law) # leaving large error bars to mimic a "rough" result (errors not displayed in plot) 
    ls_alpha2 = LeastSquares(tlist, alpha2_list, 0*tlist + 0.01, power_law) # leaving large error bars to mimic a "rough" result (errors not displayed in plot) 
    ls_alpha3 = LeastSquares(tlist, alpha3_list, 0*tlist + 0.01, power_law) # leaving large error bars to mimic a "rough" result (errors not displayed in plot) 

    m_alpha1 = Minuit(ls_alpha1, beta=1, gamma=2.5)
    m_alpha2 = Minuit(ls_alpha2, beta=1, gamma=2.5)
    m_alpha3 = Minuit(ls_alpha3, beta=1, gamma=2.5)

    m_alpha1.migrad()
    m_alpha2.migrad()
    m_alpha3.migrad()

    m_alpha1.hesse()
    m_alpha2.hesse()
    m_alpha3.hesse()

    log_tlist = np.log(tlist)

    ax.scatter(log_tlist, np.log(alpha1_list), color="#ec1010", marker='o', label=r"$\alpha_1(t, g_{12}, g_{13}, g_{23})$")
    ax.scatter(log_tlist, np.log(alpha2_list), color="#ffe600", marker='o', label=r"$\alpha_2(t, g_{12}, g_{13}, g_{23})$")
    ax.scatter(log_tlist, np.log(alpha3_list), color="#970fa3", marker='o', label=r"$\alpha_3(t, g_{12}, g_{13}, g_{23})$")
    ax.plot(log_tlist, np.log(normA_Trotter1), color="#f16608ff", label="Simulation")
    ax.plot(log_tlist, np.log(power_law(tlist, m_alpha1.values[0], m_alpha1.values[1])), color="#ec1010ab", linestyle="--")
    ax.plot(log_tlist, np.log(power_law(tlist, m_alpha2.values[0], m_alpha2.values[1])), color="#ffe60086", linestyle="--")
    ax.plot(log_tlist, np.log(power_law(tlist, m_alpha3.values[0], m_alpha3.values[1])), color="#970fa375", linestyle="--")
    ax.set_title(r"Evolution of first order additive Trotter error for $N=3$.")
    ax.set_xlabel(r"$\text{log}(t)$")
    ax.set_ylabel(r"$\text{log}\left(\left\|\mathcal{A}(t)\right\|\right)$")
    ax.set_xlim(np.min(np.log(tlist)), np.max(np.log(tlist)))

    legend = ax.legend(loc=2, frameon=True, borderaxespad=0.8, fontsize=8)
    legend.get_frame().set_facecolor('white')
    legend.get_frame().set_alpha(1.0)
    legend.get_frame().set_boxstyle("Square")

    textstr = '\n'.join([fr"$\gamma_1$: {m_alpha1.values[1]:.5f}", fr"$\gamma_2$: {m_alpha2.values[1]:.5f}", fr"$\gamma_3$: {m_alpha3.values[1]:.5f}"])   
    anchored_text = AnchoredText(textstr, loc="lower right", prop=dict(size=8), frameon=True)
    ax.add_artist(anchored_text)

    plt.grid()
    plt.tight_layout()
    fig.savefig(f'..\\plots\\3q_Trotter1_compare_bounds\\additive_err_compare.pdf', dpi=300)
    plt.close(fig)

    # Compare scaling of additive Trotter error with different choices of r at fixed evolution time T

    T = 10
    rlist = np.concatenate([np.array([1, 5]), np.linspace(10, 1000, 10)]) # ordered list

    fig, ax = plt.subplots(1, 1)
    fig.set_size_inches(16 / 2.54, 10 / 2.54)

    for r in rlist:
        tlist = np.linspace(0.01, T, int(r))
        normA_Trotter1 = np.array([np.max(np.linalg.svd(((-1j*(H_Heisenberg + E_Trotter1*t)*t).expm() - (-1j*H_Heisenberg*t).expm()).full(), compute_uv=False)) for t in tlist])

        x = (r - rlist[0]) / (rlist[-1] - rlist[0])
        x = x**4
        color = x * np.array([1,0,0]) + (1 - x) * np.array([0,0,1])

        ax.plot(tlist, normA_Trotter1, color=(
            max(0, min(1, float(color[0]))), 
            max(0, min(1, float(color[1]))), 
            max(0, min(1, float(color[2]))),
            0.5), label=f"r={r} simulation")

    #ax.plot(tlist, alpha1_list, color="#ec1010", marker='o', label=r"$\alpha_1(t, g_{12}, g_{13}, g_{23})$")
    #ax.plot(tlist, alpha2_list, color="#ffe600", marker='o', label=r"$\alpha_2(t, g_{12}, g_{13}, g_{23})$")
    #ax.plot(tlist, alpha3_list, color="#970fa3", marker='o', label=r"$\alpha_3(t, g_{12}, g_{13}, g_{23})$")
    ax.set_title(r"Evolution of first order additive Trotter error for $N=3$.")
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$\left\|\mathcal{A}(t)\right\|$")
    ax.set_xlim(np.min(tlist), np.max(tlist))

    legend = ax.legend(loc=2, frameon=True, borderaxespad=0.8, fontsize=8)
    legend.get_frame().set_facecolor('white')
    legend.get_frame().set_alpha(1.0)
    legend.get_frame().set_boxstyle("Square")

    plt.grid()
    plt.tight_layout()
    fig.savefig(f'..\\plots\\3q_Trotter1_compare_bounds\\additive_different_r.pdf', dpi=300)
    plt.close(fig)





if __name__ == "__main__":
    main()
