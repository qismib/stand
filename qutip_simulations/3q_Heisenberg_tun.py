# --------------------------- 3q_Heisenberg_tun.py --------------------------
# ---------------------------------------------------------------------------
# 

import qutip as qt
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
matplotlib.use('Agg')  

N = 3 # number of qubits

sx1 = qt.tensor(qt.sigmax(), qt.qeye(2), qt.qeye(2))
sy1 = qt.tensor(qt.sigmay(), qt.qeye(2), qt.qeye(2))
sz1 = qt.tensor(qt.sigmaz(), qt.qeye(2), qt.qeye(2))
sm1 = qt.tensor(qt.sigmam(), qt.qeye(2), qt.qeye(2))
sx2 = qt.tensor(qt.qeye(2), qt.sigmax(), qt.qeye(2))
sy2 = qt.tensor(qt.qeye(2), qt.sigmay(), qt.qeye(2))
sz2 = qt.tensor(qt.qeye(2), qt.sigmaz(), qt.qeye(2))
sm2 = qt.tensor(qt.qeye(2), qt.sigmam(), qt.qeye(2))
sx3 = qt.tensor(qt.qeye(2), qt.qeye(2), qt.sigmax())
sy3 = qt.tensor(qt.qeye(2), qt.qeye(2), qt.sigmay())
sz3 = qt.tensor(qt.qeye(2), qt.qeye(2), qt.sigmaz())
sm3 = qt.tensor(qt.qeye(2), qt.qeye(2), qt.sigmam())


# -|analog XY|-|Rx(-pi/2)|-|analog XY|-|Rx(pi/2)|-|Ry(-pi/2)|-|analog XY|-|Ry(pi/2)|-
# ----[tau]-------[ts]-------[tau]--------[ts]-------[ts]--------[tau]-------[ts]----
# -----------tau--------tau+ts-----2tau+ts----2tau+2ts----2tau+3ts---3tau+3ts---3tau+4ts

def H_XY(t, args):
    g12 = args['g12']
    g13 = args['g13'] 
    g23 = args['g23'] 

    H = g12*(sx1*sx2 + sy1*sy2) + g13*(sx1*sx3 + sy1*sy3) + g23*(sx2*sx3 + sy2*sy3)
    return H


def H_t(t, args): # TODO free hamiltonian?
    tau, ts = args["time_segments"] # tau: duration of analog gate, ts: duration of single qubit gate (TODO should I implement this using an external drive hamiltonian? YESS)
    
    # analog gate
    if t <= tau or tau+ts < t <= 2*tau + ts or 2*tau + 3*ts < t <= 3*tau + 3*ts:
        return H_XY(t, args)
    elif tau < t <= tau + ts:
        return 0.5*(-np.pi/2)*(sx1 + sx2 + sx3)
    elif 2*tau < t <= 2*tau + 2*ts:
        return 0.5*(np.pi/2)*(sx1 + sx2 + sx3)
    elif 2*tau + 2*ts < t <= 2*tau + 3*ts:
        return 0.5*(-np.pi/2)*(sy1 + sy2 + sy3)
    elif 3*tau + 3*ts < t <= 3*tau + 4*ts:
        return 0.5*(np.pi/2)*(sy1 + sy2 + sy3)
    
    print(f"Hamiltonian returned is 0 at t {t}")
    return 0


def main():
    t_analog = 1
    t_single = 0.1
    t_cycle = 3*t_analog + 4*t_single
    tsw = [t_analog, t_analog+t_single, 2*t_analog+t_single, 2*t_analog+2*t_single, 2*t_analog+3*t_single, 3*t_analog+3*t_single, 3*t_analog+4*t_single] # instants in which hamiltonian is changed i.e. a different operation is performed

    g12 = 0.1
    g13 = 0.5
    g23 = 0.1
    w1 = 5
    w2 = 4.95
    w3 = 4.90
    tlist = np.linspace(0, 3*t_analog + 4*t_single, 100)
    tlist = [t for t in tlist if t < t_cycle] # avoids t for which Hamiltonian is not defined



    H_evo = qt.QobjEvo(H_t, args={"w_vec": np.array([w1, w2, w3]), "g12": g12, "g13": g13, "g23": g23, "time_segments": (1, 0.1)})
    psi0 = qt.tensor(qt.basis(2,0), qt.basis(2,1), (qt.basis(2,0)+qt.basis(2,1))/np.sqrt(2)) # |0> @ |1> @ |+>

    res = qt.mesolve(H_evo, psi0, tlist, [], [sz1, sz2, sm1*sm1.dag(), sm2*sm2.dag(), sm3*sm3.dag()])

    n1 = res.expect[2]
    n2 = res.expect[3]
    n3 = res.expect[4]

    # plotting

    fig, ax = plt.subplots(1, 1)
    fig.set_size_inches(16 / 2.54, 10 / 2.54)

    ax.plot(tlist, n1, color="blue", label=r"$\langle \sigma_-^{(1)}\sigma_+^{(1)}\rangle$")
    ax.plot(tlist, n2, color="cyan", label=r"$\langle \sigma_-^{(2)}\sigma_+^{(2)}\rangle$")
    ax.plot(tlist, n3, color="#ae18db", label=r"$\langle \sigma_-^{(3)}\sigma_+^{(3)}\rangle$")
    plt.vlines(tsw, ymin=0,ymax=1,colors="#7B777C",ls="--")


    ax.set_title(r"Evolution of number operators expectation values")
    ax.set_xlabel("Time")
    ax.set_ylabel(r"$\langle \sigma_-\sigma_+\rangle$")
    ax.set_xlim(np.min(tlist), np.max(tlist))
    legend = ax.legend(loc=3, frameon=True, borderaxespad=0.8, fontsize=8)
    legend.get_frame().set_facecolor('white')
    legend.get_frame().set_alpha(1.0)
    legend.get_frame().set_boxstyle("Square")

    plt.grid()
    plt.tight_layout()
    fig.savefig(f'..\\plots\\3q_Heisenberg_tun\\num_op_evo.pdf', dpi=300)
    plt.close(fig)



if __name__ == "__main__":
    main()