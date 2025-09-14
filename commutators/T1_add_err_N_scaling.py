# ---------------------------- add_err_N_scaling.py ---------------------------------
#------------------------------------------------------------------------------------
# (builds up on T1_commutators_to_latex.py)
# Providing bounds for the additive Trotter error requires computing spectral norms
# of commutators of complicated Hamiltonian operators. 
# In the digital-analog setup I consider, the Hamiltonian is made up of 
# H_{XY}^N, H_{XZ}^N, H_{YZ}^N operators (see Section 5). 
# I'm interested in providing bounds for the scaling of the additive Trotter
# error as a function of the qubit number N. 
# Providing analytical bounds is challenging, so the aim of this script is 
# to perform a fit for the scaling by computing numerically the number of terms
# involved in the commutators for varying N. In order to provide the strictest 
# bound in each case I perform the computation with the same logic of the first
# calculation in Section 4.1.1.

# I aim to do N = 3,4,5,6 simulations for the first order Trotter error.
# Note: N=2 additive Trotter error is zero.

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

def line(x, alpha, beta):
    return alpha*x + beta



# ---------- create_coupling_dict(dim, uniform=True) ----------

def create_coupling_dict(N, uniform=True):
    """
    dim (int): number of qubits in the system.
    uniform (bool) (default: True): couplings are all equal in value.

    This function is designed to fit the all-to-all neutrino setup.
    """
    dic = {}
    
    if uniform:
        for i in range(N):
            for j in range(i+1, N):
                dic[f"g{i+1}{j+1}"] = 0.1
    else:
        # TODO: implement narrow cone of forward peaked neutrinos distributions (see https://arxiv.org/pdf/2207.03189)
        raise Exception("Can only generate uniform distribution of couplings.")
        pass

    return dic


# ---------- gen_Heisenberg_terms(N) ----------

def gen_Heisenberg_terms(N):
        """
        Generate the Heisenberg all-to-all Hamiltonian terms for N qubits organized
        in terms of the three Hamiltonian decomposition.

        Returns a tuple of three lists (H_{XY}^N, H_{XZ}^N, H_{YZ}^N) corresponding to the
        XX+YY, XX+ZZ, and YY+ZZ terms respectively for all pairs of qubits.
        """
        XY = [] 
        XZ = []   
        YZ = [] 

        for i in range(N):
            for j in range(i+1, N):
                gate_id = f"g{i+1}{j+1}" 

                xx_term = ["X" if k == i or k == j else "I" for k in range(N)]
                yy_term = ["Y" if k == i or k == j else "I" for k in range(N)]
                zz_term = ["Z" if k == i or k == j else "I" for k in range(N)]

                XY.append((0.5, gate_id, "".join(xx_term)))
                XY.append((0.5, gate_id, "".join(yy_term)))

                XZ.append((0.5, gate_id, "".join(xx_term)))
                XZ.append((0.5, gate_id, "".join(zz_term)))

                YZ.append((0.5, gate_id, "".join(yy_term)))
                YZ.append((0.5, gate_id, "".join(zz_term)))

        return XY, XZ, YZ

# ---------- group_by_pauli_string(comm) ----------

def group_by_pauli_string(comm):

    def is_coeff_equal(str1, str2):
        pair1 = str1.split(")(")
        pair2 = str2.split(")(")
        pair1 = [p.strip('(').strip(')').strip('g') for p in pair1]
        pair2 = [p.strip('(').strip(')').strip('g') for p in pair2]

        return sorted(pair1) == sorted(pair2)
        
    list_terms = {}
    for term in comm:
        if term[2] != "000":
            if term[2] not in list_terms:
                list_terms[term[2]] = [[term[0], term[1]]]
            else:
                cnt_not = 0 
                for i, elem in enumerate(list_terms[term[2]]):
                    if is_coeff_equal(term[1], elem[1]):
                        list_terms[term[2]][i][0] += elem[0] # refactor coefficient
                    else:
                        cnt_not +=1 
                if cnt_not == len(list_terms[term[2]]): list_terms[term[2]].append([term[0], term[1]])

    return list_terms

# ---------- compute_fixed_t_commsum(terms, t, N) ----------

def compute_fixed_t_commsum(terms, couplings, t, N):
    """
    terms (list of dicts): [{"XYZ": [[-0.5j, "(g12)(g23)"], ...], ...}, {"XYZ": [[0.5j, "(g12)(g13)"], ...], ...}]
    couplings (dict): {"g12": 0.1, "g13": 0.2, ...}
    t (float): time length of simulation
    order (int): order of Trotter formula

    return: value of first order Trotter error for a given set of couplings and time  
    """

    def eval_couplings(coeff_str, coupl_dict):
        split_str = coeff_str.split(")(")
        split_str = [p.strip(')').strip('(') for p in split_str]

        prod = 1
        for term in split_str:
            prod *= coupl_dict[term]
            
        return prod
    
    # err_vl = sum (g coeff)||Pauli string||(=1) * t^2/2
    err_vl = 0

    for comm_dict in terms:
        for _, terms in comm_dict.items():
            coeff = 0
            for term in terms:
                if term[1] != '': coeff += float(str(term[0]).strip('j')) * eval_couplings(term[1], couplings)

            err_vl += abs(coeff)

    err_vl = 0.5 * err_vl * t ** 2

    return err_vl

# ---------- color_gradient(start, end ,n)  ---------- 

def color_gradient(start, end, n):
    """Generate n hex colors from start to end."""
    start_rgb = np.array(mcolors.to_rgb(start))
    end_rgb = np.array(mcolors.to_rgb(end))
    return [ mcolors.to_hex(start_rgb + (end_rgb - start_rgb) * i/(n-1)) for i in range(n) ]


# ---------- MAIN ----------

def main():
    
    max_N = 8

    couplings = create_coupling_dict(N=max_N, uniform=True)

    ctl3 = PauliCommutators(N=3)
    ctl4 = PauliCommutators(N=4)
    ctl5 = PauliCommutators(N=5)
    ctl6 = PauliCommutators(N=6)
    ctl7 = PauliCommutators(N=7)
    ctl8 = PauliCommutators(N=8)

    A3, B3, C3 = gen_Heisenberg_terms(3)
    A4, B4, C4 = gen_Heisenberg_terms(4)
    A5, B5, C5 = gen_Heisenberg_terms(5)
    A6, B6, C6 = gen_Heisenberg_terms(6)
    A7, B7, C7 = gen_Heisenberg_terms(7)
    A8, B8, C8 = gen_Heisenberg_terms(8)

    # ----------------------------------------------------------------------------------------------------------
    # first order commutators
    A3B3 = ctl3.comm_lincombo(A3, B3)
    A3C3 = ctl3.comm_lincombo(A3, C3)
    B3C3 = ctl3.comm_lincombo(B3, C3)
    A4B4 = ctl4.comm_lincombo(A4, B4)
    A4C4 = ctl4.comm_lincombo(A4, C4)
    B4C4 = ctl4.comm_lincombo(B4, C4)
    A5B5 = ctl5.comm_lincombo(A5, B5)
    A5C5 = ctl5.comm_lincombo(A5, C5)
    B5C5 = ctl5.comm_lincombo(B5, C5)
    A6B6 = ctl6.comm_lincombo(A6, B6)
    A6C6 = ctl6.comm_lincombo(A6, C6)
    B6C6 = ctl6.comm_lincombo(B6, C6)
    A7B7 = ctl7.comm_lincombo(A7, B7)
    A7C7 = ctl7.comm_lincombo(A7, C7)
    B7C7 = ctl7.comm_lincombo(B7, C7)
    A8B8 = ctl8.comm_lincombo(A8, B8)
    A8C8 = ctl8.comm_lincombo(A8, C8)
    B8C8 = ctl8.comm_lincombo(B8, C8)

    terms_A3B3 = group_by_pauli_string(A3B3)
    terms_A3C3 = group_by_pauli_string(A3C3)
    terms_B3C3 = group_by_pauli_string(B3C3)
    terms_A4B4 = group_by_pauli_string(A4B4)
    terms_A4C4 = group_by_pauli_string(A4C4)
    terms_B4C4 = group_by_pauli_string(B4C4)
    terms_A5B5 = group_by_pauli_string(A5B5)
    terms_A5C5 = group_by_pauli_string(A5C5)
    terms_B5C5 = group_by_pauli_string(B5C5)
    terms_A6B6 = group_by_pauli_string(A6B6)
    terms_A6C6 = group_by_pauli_string(A6C6)
    terms_B6C6 = group_by_pauli_string(B6C6)
    terms_A7B7 = group_by_pauli_string(A7B7)
    terms_A7C7 = group_by_pauli_string(A7C7)
    terms_B7C7 = group_by_pauli_string(B7C7)
    terms_A8B8 = group_by_pauli_string(A8B8)
    terms_A8C8 = group_by_pauli_string(A8C8)
    terms_B8C8 = group_by_pauli_string(B8C8)

    # ----------------------------------------------------------------------------------------------------------
    # compute contributions to the upper bound for the first order Trotter error according to Propostion 9
    # in https://journals.aps.org/prx/pdf/10.1103/PhysRevX.11.011020 or equation (5.3) of Overleaf file 

    dt_list = np.arange(1,10) * 0.1

    err_vls = []
    for dt in dt_list:
        vls = []
        vls.append(compute_fixed_t_commsum([terms_A3B3, terms_A3C3, terms_B3C3], couplings, dt, 3))
        vls.append(compute_fixed_t_commsum([terms_A4B4, terms_A4C4, terms_B4C4], couplings, dt, 4)) 
        vls.append(compute_fixed_t_commsum([terms_A5B5, terms_A5C5, terms_B5C5], couplings, dt, 5)) 
        vls.append(compute_fixed_t_commsum([terms_A6B6, terms_A6C6, terms_B6C6], couplings, dt, 6)) 
        vls.append(compute_fixed_t_commsum([terms_A7B7, terms_A7C7, terms_B7C7], couplings, dt, 7)) 
        vls.append(compute_fixed_t_commsum([terms_A8B8, terms_A8C8, terms_B8C8], couplings, dt, 8)) 
        err_vls.append(vls)

    # ----------------------------------------------------------------------------------------------------------
    # fit of first order error data points with a line (in log-log scale)  

    N_vls = np.array([3, 4, 5, 6, 7, 8])
    log_N_vls = np.log(N_vls)
    log_err_vls = [np.log(vls) for vls in err_vls]
    m_objects = []

    for i, log_vls in enumerate(log_err_vls): 
        ls_obj = LeastSquares(log_N_vls, log_vls, 0*log_N_vls + 0.001, line)
        m_obj = Minuit(ls_obj, alpha=1, beta=1)
        m_obj.migrad()
        m_obj.hesse()
        m_objects.append(m_obj)

    # ----------------------------------------------------------------------------------------------------------
    # plot bound for the first order Trotter error as a function of N for various dt
    fig, ax = plt.subplots(1, 1)
    fig.set_size_inches(16 / 2.54, 10 / 2.54)

    colors = color_gradient("#ff0000", "#0000ff", len(dt_list))
    colors_line = color_gradient("#ff0000dd", "#0000ffc7", len(dt_list))

    for i, vls_dt in enumerate(log_err_vls):
        ax.scatter(log_N_vls, vls_dt, color=colors[i], marker='o', label=rf"$dt={dt_list[i]:.2f}$")
        ax.plot(log_N_vls, line(log_N_vls, m_objects[i].values[0], m_objects[i].values[1]), color=colors_line[i], linestyle="--")

    ax.set_title (r"Comparison of first order Trotter error for varying $dt$.")
    ax.set_xlabel(r"log number of qubits $log(N)$")
    ax.set_ylabel(r"$\log(\varepsilon_1)(dt)$")
    ax.set_xticks([np.log(n) for n in N_vls])
    #ax.set_xticklabels([str(np.log(v)) for v in N_vls], rotation='horizontal', fontsize=8)

    legend = ax.legend(loc=2, frameon=True, borderaxespad=0.8, fontsize=6)
    legend.get_frame().set_facecolor('white')
    legend.get_frame().set_alpha(1.0)
    legend.get_frame().set_boxstyle("Square")

    textstr =  "dt value  |  slope\n"
    textstr += "\n".join([fr"$dt={{{dt_list[i]:.2f}}}$: {m_obj.values[0]:.5f}" for i, m_obj in enumerate(m_objects)])   
    anchored_text = AnchoredText(textstr, loc="lower right", prop=dict(size=8), frameon=True)
    ax.add_artist(anchored_text)

    plt.grid()
    plt.tight_layout()
    fig.savefig(f'..\\plots\\T1_add_err_N_scaling\\first_order_N_scaling.pdf', dpi=300)
    plt.close(fig)

    # ----------------------------------------------------------------------------------------------------------
    # bound for the first order Trotter error as a function of dt for various N's
    dt_list = np.arange(1,20) * 0.1

    fig, ax = plt.subplots(1, 1)
    fig.set_size_inches(16 / 2.54, 10 / 2.54)

    err_vls = []
    err_vls.append([compute_fixed_t_commsum([terms_A3B3, terms_A3C3, terms_B3C3], couplings, dt, 3) for dt in dt_list])
    err_vls.append([compute_fixed_t_commsum([terms_A4B4, terms_A4C4, terms_B4C4], couplings, dt, 4) for dt in dt_list])
    err_vls.append([compute_fixed_t_commsum([terms_A5B5, terms_A5C5, terms_B5C5], couplings, dt, 5) for dt in dt_list])
    err_vls.append([compute_fixed_t_commsum([terms_A6B6, terms_A6C6, terms_B6C6], couplings, dt, 6) for dt in dt_list])
    err_vls.append([compute_fixed_t_commsum([terms_A7B7, terms_A7C7, terms_B7C7], couplings, dt, 7) for dt in dt_list])
    err_vls.append([compute_fixed_t_commsum([terms_A8B8, terms_A8C8, terms_B8C8], couplings, dt, 8) for dt in dt_list])

    for i, vls in enumerate(err_vls):
        ax.scatter(dt_list, vls, color=colors[i], marker='o', label=rf"$N={i+3}$")

    ax.set_title (r"Comparison of first order Trotter error in terms of $dt$.")
    ax.set_xlabel(r"$dt$")
    ax.set_ylabel(r"$\varepsilon_1(dt)$")

    legend = ax.legend(loc=2, frameon=True, borderaxespad=0.8, fontsize=6)
    legend.get_frame().set_facecolor('white')
    legend.get_frame().set_alpha(1.0)
    legend.get_frame().set_boxstyle("Square")

    plt.grid()
    plt.tight_layout()
    fig.savefig(f'..\\plots\\T1_add_err_N_scaling\\first_order_dt_scaling.pdf', dpi=300)
    plt.close(fig)

    # ----------------------------------------------------------------------------------------------------------
    # simulate first order evolution for the N=3 system and compare with the numerically computed bound
    dt_list = np.arange(1,20) * 0.1

    sx1 = qt.tensor(qt.sigmax(), qt.qeye(2), qt.qeye(2))
    sy1 = qt.tensor(qt.sigmay(), qt.qeye(2), qt.qeye(2))
    sz1 = qt.tensor(qt.sigmaz(), qt.qeye(2), qt.qeye(2))

    sx2 = qt.tensor(qt.qeye(2), qt.sigmax(), qt.qeye(2))
    sy2 = qt.tensor(qt.qeye(2), qt.sigmay(), qt.qeye(2))
    sz2 = qt.tensor(qt.qeye(2), qt.sigmaz(), qt.qeye(2))

    sx3 = qt.tensor(qt.qeye(2), qt.qeye(2), qt.sigmax())
    sy3 = qt.tensor(qt.qeye(2), qt.qeye(2), qt.sigmay())
    sz3 = qt.tensor(qt.qeye(2), qt.qeye(2), qt.sigmaz())

    g12 = couplings["g12"]
    g13 = couplings["g13"]
    g23 = couplings["g23"]

    H_Heisenberg = g12*(sx1*sx2 + sy1*sy2 + sz1*sz2) + g13*(sx1*sx3 + sy1*sy3 + sz1*sz3) + g23*(sx2*sx3 + sy2*sy3 + sz2*sz3)

    E_Trotter1  = 0.25 * (sx1*sz2*sy3 + sz1*sx2*sy3) * (g12*g13+g12*g23-3*g13*g23)
    E_Trotter1 += 0.25 * (sx1*sy2*sz3 + sz1*sy2*sz3) * (g12*g13+g13*g23-3*g12*g23)
    E_Trotter1 += 0.25 * (sy1*sx2*sz3 + sy1*sz2*sx3) * (g12*g23+g13*g23-3*g12*g13)

    normA_Trotter1 = np.array([np.max(np.linalg.svd(((-1j*(H_Heisenberg + E_Trotter1*t)*t).expm() - (-1j*H_Heisenberg*t).expm()).full(), compute_uv=False)) for t in dt_list])
    errvls3 = [compute_fixed_t_commsum([terms_A3B3, terms_A3C3, terms_B3C3], couplings, dt, 3) for dt in dt_list]

    fig, ax = plt.subplots(1, 1)
    fig.set_size_inches(16 / 2.54, 10 / 2.54)

    ax.scatter(dt_list, errvls3, color="#0026FF", marker='o', label=rf"Numerical bound")
    ax.scatter(dt_list, normA_Trotter1, color="#F700FF", marker='o', label=rf"Simulation")

    ax.set_title (r"Compare simulation and numerical bounds for first order Trotter error for N=3.")
    ax.set_xlabel(r"$dt$")
    ax.set_ylabel(r"$\varepsilon_1(dt)$")

    legend = ax.legend(loc=2, frameon=True, borderaxespad=0.8, fontsize=6)
    legend.get_frame().set_facecolor('white')
    legend.get_frame().set_alpha(1.0)
    legend.get_frame().set_boxstyle("Square")

    plt.grid()
    plt.tight_layout()
    fig.savefig(f'..\\plots\\T1_add_err_N_scaling\\N3_scalings.pdf', dpi=300)
    plt.close(fig)

if __name__ == "__main__":
    main()