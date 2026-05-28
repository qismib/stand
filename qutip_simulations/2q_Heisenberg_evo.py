# ------------------------- 2q_Heisenberg_evo.py ----------------------------
# ---------------------------------------------------------------------------
# Analysis of a 2 qubit system evolving under a time dependent version of the
# Heisenberg hamiltonian.

import qutip as qt
import numpy as np
import sys; sys.path.append("../classes")
from PlotterTool import BasicPlotter

import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg')  



def Rx(theta):
    return (-1j * theta/2 * qt.sigmax()).expm()

def Ry(theta):
    return (-1j * theta/2 * qt.sigmay()).expm()

def Rz(theta):
    return (-1j * theta/2 * qt.sigmaz()).expm()

def H_t(t, args):
        # Heisenberg
        g = args['g']
        w = args["w"]
        H0 = 0.5*w * (qt.tensor(qt.sigmaz(), qt.qeye(2)) + qt.tensor(qt.qeye(2), qt.sigmaz()))
        H = H0 + g*np.abs((np.cos(0.1*w*t)))*(qt.tensor(qt.sigmax(), qt.sigmax()) + qt.tensor(qt.sigmay(), qt.sigmay()) + qt.tensor(qt.sigmaz(), qt.sigmaz()))
        return H


# ------------------- MAIN -------------------

def main():
    g = 0.1
    w = 5 # qubits are at resonance
    tlist = np.linspace(-20, 20, 1000)

    sz1 = qt.tensor(qt.sigmaz(), qt.qeye(2))
    sm1 = qt.tensor(qt.sigmam(), qt.qeye(2))
    sz2 = qt.tensor(qt.qeye(2), qt.sigmaz())
    sm2 = qt.tensor(qt.qeye(2), qt.sigmam())

    H_evo = qt.QobjEvo(H_t, args={"w": w, "g": g})
    psi0 = qt.tensor(qt.basis(2,0), qt.basis(2,1)) # |0> @ |1>

    res = qt.mesolve(H_evo, psi0, tlist, [], [sz1, sz2, sm1*sm1.dag(), sm2*sm2.dag()])

    plotter = BasicPlotter(tlist, res.expect[0:2], ("Time", ""))
    plotter.labels = ["$\\langle \\sigma_z^1 \\rangle$", "$\\langle \\sigma_z^2 \\rangle$"]
    plotter.title = "$\\langle\\sigma_z\\rangle$ with $g=g(t)$."
    plotter.file_name = "../plots/2q_Heisenberg_evo/expval_sz_gt.pdf"
    plotter.save_plot()

    plotter = BasicPlotter(tlist, res.expect[2:], ("Time", ""))
    plotter.labels = ["$\\langle \\sigma_-^1\\sigma_+^1 \\rangle$", "$\\langle \\sigma_-^2\\sigma_+^2 \\rangle$"]
    plotter.title = "$\\langle\\sigma_-\\sigma_+\\rangle$ with $g=g(t)$."
    plotter.file_name = "../plots/2q_Heisenberg_evo/expval_numop_gt.pdf"
    plotter.save_plot()

    




if __name__ == "__main__":
    main()

