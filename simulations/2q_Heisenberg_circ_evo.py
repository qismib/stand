import qutip as qt
import numpy as np
import sys; sys.path.append("../classes")

import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg')  



def Rx(theta):
    return (-1j * theta/2 * qt.sigmax()).expm()

def Ry(theta):
    return (-1j * theta/2 * qt.sigmay()).expm()

def Rz(theta):
    return (-1j * theta/2 * qt.sigmaz()).expm()



# ------------------- MAIN -------------------

def main():
    n_step = 100 # number of repeated applications of U_XY
    g = 0.1
    tau = 250 # 250ns

    sx1 = qt.tensor(qt.sigmax(), qt.qeye(2))
    sy1 = qt.tensor(qt.sigmay(), qt.qeye(2))
    sz1 = qt.tensor(qt.sigmaz(), qt.qeye(2))
    sx2 = qt.tensor(qt.qeye(2), qt.sigmax())
    sy2 = qt.tensor(qt.qeye(2), qt.sigmay())
    sz2 = qt.tensor(qt.qeye(2), qt.sigmaz())

    H_XY = 0.5*g*(sx1*sx2 + sy1*sy2)

    U_XY = (-1j*H_XY*tau).expm()
    U_XZ = qt.tensor(Rx(np.pi/2), Rx(np.pi/2)) * U_XY * qt.tensor(Rx(-np.pi/2), Rx(-np.pi/2))
    U_YZ = qt.tensor(Ry(np.pi/2), Ry(np.pi/2)) * U_XY * qt.tensor(Ry(-np.pi/2), Ry(-np.pi/2))

    U_full = U_YZ * U_XZ * U_XY

    psi0 = qt.tensor(qt.basis(2,1), qt.basis(2,0))

    sz1_values = np.zeros(n_step, dtype=complex)
    sz2_values = np.zeros(n_step, dtype=complex)

    sz1_values[0] = psi0.dag() * sz1 * psi0
    sz2_values[0] = psi0.dag() * sz2 * psi0

    for i in range(1, n_step):
        psi_i = U_full * psi0
        sz1_values[i] = psi_i.dag() * sz1 * psi_i
        sz2_values[i] = psi_i.dag() * sz2 * psi_i
        psi0 = psi_i

    plt.plot(np.linspace(0, n_step*tau, n_step), sz1_values, label="$\\langle \\sigma^{x_1} \\rangle$")
    plt.plot(np.linspace(0, n_step*tau, n_step), sz2_values, label="$\\langle \\sigma^{x_2} \\rangle$")
    plt.xlabel("Time")
    plt.ylabel("$\\langle \\sigma_x^{1,2} \\rangle$")

    plt.grid()
    plt.legend(loc=1)

    plt.savefig("../plots/2q_Heisenberg_cir_evo/test.pdf")


    



if __name__ == "__main__":
    main()






