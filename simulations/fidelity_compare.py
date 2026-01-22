# ------------------------ fidelity_compare.py -----------------------
#---------------------------------------------------------------------
# Compare fidelity values of single-step algorithm as we let the 
# parameters of the system vary. Compare the fully digital and the
# digital-analog approach.
# My reference for the digital approach is:
# https://journals.aps.org/prd/pdf/10.1103/PhysRevD.107.023007
# I compare the first order versions of the algorithms.

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors  
matplotlib.use('Agg')  

def fidelityDigital(N, args):
    return (args["F_1q"] ** args["N1D"](N)) * (args["F_Uq"] ** args["NUq"](N)) * (args["F_zz"] ** args["Nzz"](N))

def fidelityDAnalog(N, args):
    return (args["F_1q"] ** args["N1DA"](N)) * (args["F_xy"](N) ** args["Nxy"])

def color_gradient(start, end, n):
    """Generate n hex colors from start to end."""
    start_rgb = np.array(mcolors.to_rgb(start))
    end_rgb = np.array(mcolors.to_rgb(end))
    return [mcolors.to_hex(start_rgb + (end_rgb - start_rgb) * i/(n-1)) for i in range(n)]

def fidelityEvaluator(values_N, args):
    fid_dig = np.zeros(len(values_N))
    fid_dan = np.zeros(len(values_N))
    fid_ratio = np.zeros(len(values_N))

    for N in values_N:
        fid_dig[N-values_N[0]] = fidelityDigital(N, args)
        fid_dan[N-values_N[0]] = fidelityDAnalog(N, args)
        fid_ratio[N-values_N[0]] = fid_dan[N-values_N[0]] / fid_dig[N-values_N[0]] # NOTE: I look for F_DA / F_D > 1 s.t. F_DA > F_D

    return fid_dig, fid_dan, fid_ratio


# def main():
#     params = {
#         # fidelities
#         "F_1q": 1-1e-3, # fidelity of 1q rotations gates
#         "F_Uq": 1-1e-3, # fidelity of Uq gates
#         "F_zz": 1-1e-2,  # fidelity of 2 qubit ZZ gates
#         "F_xy": lambda N: 0.99 - 0.025*(N-1), # fidelity of analog XY gate
# 
#         # number of gates per single Trotter step
#         "N1D" : lambda N: 3*N*(N-1)/2 + N,
#         "N1DA": lambda N: 4*N + N,       
#         "NUq" : lambda N: 6*N*(N-1),
#         "Nzz" : lambda N: 3*N*(N-1)/2,
#         "Nxy" : 3
#     }
# 
#     N_max = 10
#     values_N = np.arange(2, N_max+1)
# 
#     fig, ax = plt.subplots()
#     colors = color_gradient("#ff0000", "#0000ff", len(values_N))
# 
#     for N in values_N:
#         ax.scatter(N, fidelityDAnalog(N, params)/fidelityDigital(N, params), color=colors[N-2], label=f"N = {N}")
# 
#     ax.set_title("Compare fidelities of digital and digital-analog approach.")
# 
#     ax.set_xlabel("N")
#     ax.set_ylabel(r"$\mathcal{F}_\text{DA}/\mathcal{F}_\text{digital} (N)$")
#     legend = ax.legend(loc=2, frameon=True, borderaxespad=0.8, fontsize=6)
#     legend.get_frame().set_facecolor('white')
#     legend.get_frame().set_alpha(1.0)
#     legend.get_frame().set_boxstyle("Square")
#     plt.grid()
#     plt.tight_layout()
#     plt.savefig('..\\plots\\fidelity_compare\\T1_ratio_fidelities_digital_and_digitalanalog.pdf', dpi=300)
#     plt.close(fig)
# 
#     fig, ax = plt.subplots()
#     colors = color_gradient("#ff0000", "#0000ff", N_max-1)
# 
#     fid_digital = np.zeros(len(values_N))
#     fid_danalog = np.zeros(len(values_N))
# 
#     for N in values_N:
#         fid_digital[N-2] = fidelityDigital(N, params)
#         fid_danalog[N-2] = fidelityDAnalog(N, params)
# 
#     ax.plot(values_N, fid_danalog, label="DAnalog")
#     ax.plot(values_N, fid_digital, label="Digital")
# 
#     ax.set_title("Compare fidelities of digital and digital-analog approach.")
#     ax.set_xlabel("N")
#     ax.set_ylabel(r"$\mathcal{F}(N)$")
#     legend = ax.legend(loc=2, frameon=True, borderaxespad=0.8, fontsize=6)
#     legend.get_frame().set_facecolor('white')
#     legend.get_frame().set_alpha(1.0)
#     legend.get_frame().set_boxstyle("Square")
#     plt.grid()
#     plt.tight_layout()
#     plt.savefig('..\\plots\\fidelity_compare\\T1_plot_fidelities.pdf', dpi=300)
#     plt.close(fig)

def main():
    penalty_values = np.arange(10, 50) * 1e-3 

    N_max = 10
    values_N = np.arange(2, N_max+1)

    if 0 in values_N or N_max <= 0: raise Exception("The choose values of N are not valid.")

    values_N_advantage = []

    for penalty in penalty_values:

        params = {
            # fidelities
            "F_1q": 1-1e-3, # fidelity of 1q rotations gates
            "F_Uq": 1-1e-3, # fidelity of Uq gates
            "F_zz": 1-1e-2,  # fidelity of 2 qubit ZZ gates
            "F_xy": lambda N: 0.99 - penalty*(N-1), # fidelity of analog XY gate

            # number of gates per single Trotter step
            "N1D" : lambda N: 3*N*(N-1)/2 + N,
            "N1DA": lambda N: 4*N + N,       
            "NUq" : lambda N: 6*N*(N-1),
            "Nzz" : lambda N: 3*N*(N-1)/2,
            "Nxy" : 3
        }

        _, _, fidelityRatio = fidelityEvaluator(values_N, params)

        for i, value in enumerate(fidelityRatio): 
            if value >= 1: 
                values_N_advantage.append(i + values_N[0])
                break



    fig, ax = plt.subplots()
    colors = color_gradient("#ff0000", "#0000ff", len(penalty_values))      

    for i, penalty in enumerate(penalty_values):
        ax.scatter(penalty, values_N_advantage[i], color=colors[i])

    ax.set_title("Advantage DA vs. D as a function of penalty on analog block.")
    ax.set_xlabel("penalty")
    ax.set_ylabel(r"$N$ value")
    legend = ax.legend(loc=2, frameon=True, borderaxespad=0.8, fontsize=6)
    legend.get_frame().set_facecolor('white')
    legend.get_frame().set_alpha(1.0)
    legend.get_frame().set_boxstyle("Square")
    plt.grid()
    plt.tight_layout()
    plt.savefig('..\\plots\\fidelity_compare\\T1_Navdvantage_vs_penalty_analog_gate.pdf', dpi=300)
    plt.close(fig)



if __name__ == "__main__":
    main()