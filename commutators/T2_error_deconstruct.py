# ---------------------------- T2error_deconstruct.py ------------------------------
#------------------------------------------------------------------------------------


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

def group_by_pauli_string(comm, N):

    def is_coeff_equal(str1, str2):
        pair1 = str1.split(")(")        
        pair2 = str2.split(")(")
        pair1 = [p.strip('(').strip(')').strip('g') for p in pair1]
        pair2 = [p.strip('(').strip(')').strip('g') for p in pair2]

        return sorted(pair1) == sorted(pair2)
        
    list_terms = {}
    for term in comm:
        if term[2] != "0"*N:
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

def compute_fixed_t_commsum_T2(onetwo, twofour, couplings, t, N):
    """
    onetwo  (list of dicts): [{"XYZ": [[-0.5j, "(g12)(g23)"], ...], ...}, {"XYZ": [[0.5j, "(g12)(g13)"], ...], ...}]; these are the terms that contribute with 1/12 coeff
    twofour (list of dicts): [{"XYZ": [[-0.5j, "(g12)(g23)"], ...], ...}, {"XYZ": [[0.5j, "(g12)(g13)"], ...], ...}]; these are the terms that contribute with 1/24 coeff
    couplings (dict): {"g12": 0.1, "g13": 0.2, ...}
    t (float): time length of simulation
    order (int): order of Trotter formula

    Implements second order Trotter formula for the Trotter error.

    return: value of first order Trotter error for a given set of couplings and time.  
    """
    def format_couplings(s):
        coupl = s.replace('(', '').replace(')','')
        coupl = coupl.split('g')
        coupl = ['g'+pair for pair in coupl if pair != '']
        return coupl


    def eval_couplings(coeff_str, coupl_dict):
        split_str = format_couplings(coeff_str)

        prod = 1
        for term in split_str:
            prod *= coupl_dict[term]
            
        return prod
    
    # err_vl = 1/12 * sum (g coeff)||Pauli string||(=1) * dt^3 + 1/24 * sum (g coeff)||Pauli string||(=1) * dt^3 (TODO check)
    err_vl = 0

    err_12 = 0
    for comm_dict in onetwo:
        for _, terms in comm_dict.items():
            coeff = 0
            for term in terms:
                if term[1] != '': coeff += term[0].real * eval_couplings(term[1], couplings)
            err_12 += abs(coeff)

    err_vl += err_12 * t**3 / 12

    err_24 = 0
    for comm_dict in twofour:
        for _, terms in comm_dict.items():
            coeff = 0
            for term in terms:
                if term[1] != '': coeff += term[0].real * eval_couplings(term[1], couplings)
            err_24 += abs(coeff)

    err_vl += err_24 * t**3 / 24

    return err_vl

# ---------- color_gradient(start, end, n)  ---------- 

def color_gradient(start, end, n):
    """
    start (hex str): starting color value
    end   (hex str): ending color value
    n (int): steps

    Generate n hex colors from start to end.
    """
    start_rgb = np.array(mcolors.to_rgb(start))
    end_rgb = np.array(mcolors.to_rgb(end))
    return [mcolors.to_hex(start_rgb + (end_rgb - start_rgb) * i/(n-1)) for i in range(n)]


# ---------- MAIN ----------

def main():
    
    N_max = 6

    couplings = create_coupling_dict(N=N_max, uniform=True)

    ctl3 = PauliCommutators(N=3)
    ctl4 = PauliCommutators(N=4)
    ctl5 = PauliCommutators(N=5)
    ctl6 = PauliCommutators(N=6)

    A3, B3, C3 = gen_Heisenberg_terms(3)
    A4, B4, C4 = gen_Heisenberg_terms(4)
    A5, B5, C5 = gen_Heisenberg_terms(5)
    A6, B6, C6 = gen_Heisenberg_terms(6)

    # ----------------------------------------------------------------------------------------------------------
    # first order commutators
    # NOTE: could do without redefining since - sign difference does not affect the result
    #       (absolute value is taken at the end)
    A3B3 = ctl3.comm_lincombo(A3, B3)
    A3C3 = ctl3.comm_lincombo(A3, C3)
    B3C3 = ctl3.comm_lincombo(B3, C3)
    B3A3 = ctl3.comm_lincombo(B3, A3)
    C3A3 = ctl3.comm_lincombo(C3, A3)
    C3B3 = ctl3.comm_lincombo(C3, B3)

    A4B4 = ctl4.comm_lincombo(A4, B4)
    A4C4 = ctl4.comm_lincombo(A4, C4)
    B4C4 = ctl4.comm_lincombo(B4, C4)
    B4A4 = ctl4.comm_lincombo(B4, A4)
    C4A4 = ctl4.comm_lincombo(C4, A4)
    C4B4 = ctl4.comm_lincombo(C4, B4)

    A5B5 = ctl5.comm_lincombo(A5, B5)
    A5C5 = ctl5.comm_lincombo(A5, C5)
    B5C5 = ctl5.comm_lincombo(B5, C5)
    B5A5 = ctl5.comm_lincombo(B5, A5)
    C5A5 = ctl5.comm_lincombo(C5, A5)
    C5B5 = ctl5.comm_lincombo(C5, B5)

    A6B6 = ctl6.comm_lincombo(A6, B6)
    A6C6 = ctl6.comm_lincombo(A6, C6)
    B6C6 = ctl6.comm_lincombo(B6, C6)
    B6A6 = ctl6.comm_lincombo(B6, A6)
    C6A6 = ctl6.comm_lincombo(C6, A6)
    C6B6 = ctl6.comm_lincombo(C6, B6)   

    # 1/12
    B3_B3A3 = ctl3.comm_lincombo(B3, B3A3)
    B3_C3A3 = ctl3.comm_lincombo(B3, C3A3)
    C3_B3A3 = ctl3.comm_lincombo(C3, B3A3)
    C3_C3A3 = ctl3.comm_lincombo(C3, C3A3)
    C3_C3B3 = ctl3.comm_lincombo(C3, C3B3)
    # 1/24
    A3_A3B3 = ctl3.comm_lincombo(A3, A3B3)
    A3_A3C3 = ctl3.comm_lincombo(A3, A3C3)
    B3_B3C3 = ctl3.comm_lincombo(B3, B3C3)
    
    # 1/12
    B4_B4A4 = ctl4.comm_lincombo(B4, B4A4)
    B4_C4A4 = ctl4.comm_lincombo(B4, C4A4)
    C4_B4A4 = ctl4.comm_lincombo(C4, B4A4)
    C4_C4A4 = ctl4.comm_lincombo(C4, C4A4)
    C4_C4B4 = ctl4.comm_lincombo(C4, C4B4)
    # 1/24
    A4_A4B4 = ctl4.comm_lincombo(A4, A4B4)
    A4_A4C4 = ctl4.comm_lincombo(A4, A4C4)
    B4_B4C4 = ctl4.comm_lincombo(B4, B4C4)

    # 1/12
    B5_B5A5 = ctl5.comm_lincombo(B5, B5A5)
    B5_C5A5 = ctl5.comm_lincombo(B5, C5A5)
    C5_B5A5 = ctl5.comm_lincombo(C5, B5A5)
    C5_C5A5 = ctl5.comm_lincombo(C5, C5A5)
    C5_C5B5 = ctl5.comm_lincombo(C5, C5B5)
    # 1/24
    A5_A5B5 = ctl5.comm_lincombo(A5, A5B5)
    A5_A5C5 = ctl5.comm_lincombo(A5, A5C5)
    B5_B5C5 = ctl5.comm_lincombo(B5, B5C5)

    # 1/12
    B6_B6A6 = ctl6.comm_lincombo(B6, B6A6)
    B6_C6A6 = ctl6.comm_lincombo(B6, C6A6)
    C6_B6A6 = ctl6.comm_lincombo(C6, B6A6)
    C6_C6A6 = ctl6.comm_lincombo(C6, C6A6)
    C6_C6B6 = ctl6.comm_lincombo(C6, C6B6)
    # 1/24
    A6_A6B6 = ctl6.comm_lincombo(A6, A6B6)
    A6_A6C6 = ctl6.comm_lincombo(A6, A6C6)
    B6_B6C6 = ctl6.comm_lincombo(B6, B6C6)

    # N=3
    terms_B3_C3A3 = group_by_pauli_string(B3_C3A3, N=3)
    terms_B3_B3A3 = group_by_pauli_string(B3_B3A3, N=3)
    terms_C3_C3A3 = group_by_pauli_string(C3_C3A3, N=3)
    terms_C3_B3A3 = group_by_pauli_string(C3_B3A3, N=3)
    terms_A3_A3B3 = group_by_pauli_string(A3_A3B3, N=3)
    terms_C3_C3B3 = group_by_pauli_string(C3_C3B3, N=3)
    terms_B3_B3C3 = group_by_pauli_string(B3_B3C3, N=3)
    terms_A3_A3C3 = group_by_pauli_string(A3_A3C3, N=3)

    # N=4
    terms_B4_C4A4 = group_by_pauli_string(B4_C4A4, N=4)
    terms_B4_B4A4 = group_by_pauli_string(B4_B4A4, N=4)
    terms_C4_C4A4 = group_by_pauli_string(C4_C4A4, N=4)
    terms_C4_B4A4 = group_by_pauli_string(C4_B4A4, N=4)
    terms_A4_A4B4 = group_by_pauli_string(A4_A4B4, N=4)
    terms_C4_C4B4 = group_by_pauli_string(C4_C4B4, N=4)
    terms_B4_B4C4 = group_by_pauli_string(B4_B4C4, N=4)
    terms_A4_A4C4 = group_by_pauli_string(A4_A4C4, N=4)
    
    # N=5
    terms_B5_C5A5 = group_by_pauli_string(B5_C5A5, N=5)
    terms_B5_B5A5 = group_by_pauli_string(B5_B5A5, N=5)
    terms_C5_C5A5 = group_by_pauli_string(C5_C5A5, N=5)
    terms_C5_B5A5 = group_by_pauli_string(C5_B5A5, N=5)
    terms_A5_A5B5 = group_by_pauli_string(A5_A5B5, N=5)
    terms_C5_C5B5 = group_by_pauli_string(C5_C5B5, N=5)
    terms_B5_B5C5 = group_by_pauli_string(B5_B5C5, N=5)
    terms_A5_A5C5 = group_by_pauli_string(A5_A5C5, N=5)

    # N=6
    terms_B6_C6A6 = group_by_pauli_string(B6_C6A6, N=6)
    terms_B6_B6A6 = group_by_pauli_string(B6_B6A6, N=6)
    terms_C6_C6A6 = group_by_pauli_string(C6_C6A6, N=6)
    terms_C6_B6A6 = group_by_pauli_string(C6_B6A6, N=6)
    terms_A6_A6B6 = group_by_pauli_string(A6_A6B6, N=6)
    terms_C6_C6B6 = group_by_pauli_string(C6_C6B6, N=6)
    terms_B6_B6C6 = group_by_pauli_string(B6_B6C6, N=6)
    terms_A6_A6C6 = group_by_pauli_string(A6_A6C6, N=6)

    # ----------------------------------------------------------------------------------------------------------
    # compute contributions to the upper bound for the second order Trotter error according to Propostion 9
    # in https://journals.aps.org/prx/pdf/10.1103/PhysRevX.11.011020 or equation (5.3) of Overleaf file 

    dt_list = np.arange(1,10) * 0.1

    terms12_3 = [terms_B3_B3A3, terms_B3_C3A3, terms_C3_B3A3, terms_C3_C3A3, terms_C3_C3B3]
    terms12_4 = [terms_B4_B4A4, terms_B4_C4A4, terms_C4_B4A4, terms_C4_C4A4, terms_C4_C4B4] 
    terms12_5 = [terms_B5_B5A5, terms_B5_C5A5, terms_C5_B5A5, terms_C5_C5A5, terms_C5_C5B5] 
    terms12_6 = [terms_B6_B6A6, terms_B6_C6A6, terms_C6_B6A6, terms_C6_C6A6, terms_C6_C6B6] 

    terms24_3 = [terms_A3_A3B3, terms_A3_A3C3, terms_B3_B3C3]
    terms24_4 = [terms_A4_A4B4, terms_A4_A4C4, terms_B4_B4C4]
    terms24_5 = [terms_A5_A5B5, terms_A5_A5C5, terms_B5_B5C5]
    terms24_6 = [terms_A6_A6B6, terms_A6_A6C6, terms_B6_B6C6]

    # ----------------------------------------------------------------------------------------------------------
    # Estimate Trotter error vs. r dependance
    # Strategy: evaluate the total Trotter error after a simulation of total duration T for simulations 
    # performed using different numbers of Trotter steps to simulate th e dynamics.

    N = 3
    T = 5
    r_vls = np.arange(1, 21) * 5

    fig, ax = plt.subplots(1, 1)
    fig.set_size_inches(16 / 2.54, 10 / 2.54)

    colors = color_gradient("#ff0000", "#0000ff", len(r_vls))

    for i, r in enumerate(r_vls):
        dt_list = np.linspace(0, T, r)
        err_vls_r = []
        for dt in dt_list: err_vls_r.append(compute_fixed_t_commsum_T2(terms12_3, terms24_3, couplings, dt, N))

        err_vls_r_scaled = [np.log(x) for x in err_vls_r if x!=0]
        dt_list_log = [np.log(x) for x in dt_list if x!=0]
        ax.scatter(dt_list_log, err_vls_r_scaled, color=colors[i], marker='o', label=f"r={r}")

        ax.set_title (r"Comparison of seoond order Trotter error for varying $r$.")
        ax.set_xlabel(r"$log(r)$")
        ax.set_ylabel(r"$\log(\varepsilon_2)(dt)$")

        legend = ax.legend(loc=2, frameon=True, borderaxespad=0.8, fontsize=6)
        legend.get_frame().set_facecolor('white')
        legend.get_frame().set_alpha(1.0)
        legend.get_frame().set_boxstyle("Square")

    plt.grid()
    plt.tight_layout()
    fig.savefig(f'..\\plots\\T2_error_deconstruct\\second_order_scaling.pdf', dpi=300)
    plt.close(fig)

    # ----------------------------------------------------------------------------------------------------------
    T = 5.0
    N = 3
    r_vls = np.array([2,4,8,16,32,64,128])
    err_vs_r = []

    for r in r_vls:
        dt = T / float(r)               # step size used in the trotter simulation
        val = compute_fixed_t_commsum_T2(terms12_3, terms24_3, couplings, dt, N)
        err_vs_r.append(val)

    err_vs_r = np.array(err_vs_r)

    # Remove zero or negative values (if any) before taking logs
    mask = err_vs_r > 0
    r_fit = r_vls[mask]
    err_fit = err_vs_r[mask]

    # Fit log(err) = A + B * log(r) using LeastSquares + iminuit
    lx = np.log(r_fit)
    ly = np.log(err_fit)

    # Use your line(x, alpha, beta) where line returns alpha*x + beta (x is log r)
    # define LeastSquares cost for iminuit (x,y are globals here)
    cost = LeastSquares(lx, ly, ly*0 + 0.001, line)   # third arg is yerr (0 → unweighted)
    m = Minuit(cost, alpha=-2.0, beta=0.0)  # initial guess: slope ~ -2
    m.migrad()
    m.hesse()

    alpha_hat = m.values['alpha']
    beta_hat  = m.values['beta']       # intercept in log space
    alpha_err = m.errors['alpha']
    C_hat = np.exp(beta_hat)

    print(f"fit: log(err) = {alpha_hat:.3f} * log(r) + {beta_hat:.3f}")
    print(f"so err ~ C r^{alpha_hat:.3f} with C={C_hat:.3e} and σ(alpha)={alpha_err:.3f}")

    # Make a log-log plot with fit
    fig, ax = plt.subplots()
    ax.scatter(r_fit, err_fit, label='data')
    r_plot = np.logspace(np.log10(r_fit.min()), np.log10(r_fit.max()), 200)
    ax.plot(r_plot, np.exp(line(np.log(r_plot), alpha_hat, beta_hat)), 
            label=f'fit slope={alpha_hat:.3f}', linestyle='--')
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlabel('r (Trotter steps)')
    ax.set_ylabel('Trotter error (estimate)')
    ax.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig('..\\plots\\T2_error_deconstruct\\error_vs_r_corrected.pdf', dpi=300)
    plt.close(fig)




    #err_vls = []
    #for dt in dt_list:
    #    vls = []
    #    vls.append(compute_fixed_t_commsum_T2(terms12_3, terms24_3, couplings, dt, 3))
    #    vls.append(compute_fixed_t_commsum_T2(terms12_4, terms24_4, couplings, dt, 4)) 
    #    vls.append(compute_fixed_t_commsum_T2(terms12_5, terms24_5, couplings, dt, 5)) 
    #    vls.append(compute_fixed_t_commsum_T2(terms12_6, terms24_6, couplings, dt, 6)) 
    #    err_vls.append(vls)

    #fig, ax = plt.subplots(1, 1)
    #fig.set_size_inches(16 / 2.54, 10 / 2.54)

    #for i, elem in enumerate(err_vls):
    #    ax.scatter(dt_list, elem, color="", marker='o')
    #    #ax.plot(log_N_vls, line(log_N_vls, m_objects[i].values[0], m_objects[i].values[1]), color=colors_line[i], linestyle="--")

    #ax.set_title (r"Comparison of seoond order Trotter error for varying $dt$.")
    #ax.set_xlabel(r"log number of qubits $log(N)$")
    #ax.set_ylabel(r"$\log(\varepsilon_1)(dt)$")

    #legend = ax.legend(loc=2, frameon=True, borderaxespad=0.8, fontsize=6)
    #legend.get_frame().set_facecolor('white')
    #legend.get_frame().set_alpha(1.0)
    #legend.get_frame().set_boxstyle("Square")

    #plt.grid()
    #plt.tight_layout()
    #fig.savefig(f'..\\plots\\T2_error_deconstruct\\second_order_N_scaling.pdf', dpi=300)
    #plt.close(fig)


    

if __name__ == "__main__":
    main()