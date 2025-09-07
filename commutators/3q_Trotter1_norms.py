import numpy as np
import sys; sys.path.append("../classes")
from CommutatorTool import CommutatorNorms




# ------------------- MAIN -------------------

def main():
    # Calculations relevant to paragraph 4.1.1.
    # I separately compute the spectral norm of each sum of 6 Pauli strings.
    N = 3 # number of qubits

    tool_3q = CommutatorNorms(N)

    # strings of three distinct Pauli matrices
    XYZ = tool_3q.pauli_tensor("XYZ")
    XZY = tool_3q.pauli_tensor("XZY")
    YXZ = tool_3q.pauli_tensor("YXZ")
    YZX = tool_3q.pauli_tensor("YZX")
    ZXY = tool_3q.pauli_tensor("ZXY")
    ZYX = tool_3q.pauli_tensor("ZYX")
    paulis = [XYZ, XZY, YXZ, YZX, ZXY, ZYX]

    # lists are provided in the same order of definition of strings of 3 distinct Pauli (True -> +, False -> -)
    # TODO compute this automatic
    list1 = [True, True, False, False, False, False]
    list2 = [True, False, True, False, True, True]
    list3 = [False, False, False, True, False, True]
    list4 = [True, True, False, False, True, True]
    list5 = [False, True, True, True, True, False]
    list6 = [True, False, True, True, False, True]
    list7 = [False, False, False, False, True, True]
    list8 = [False, True, False, True, False, False]
    list9 = [True, False, True, False, False, False]

    def sum_with_sign(list):
        sum = tool_3q.zero_obj()
        for i, elem in enumerate(list):
            if elem: 
                sum += paulis[i]
            else:
                sum -= paulis[i]
        return sum

    sum1 = sum_with_sign(list1)
    sum2 = sum_with_sign(list2)
    sum3 = sum_with_sign(list3)
    sum4 = sum_with_sign(list4)
    sum5 = sum_with_sign(list5)
    sum6 = sum_with_sign(list6)
    sum7 = sum_with_sign(list7)
    sum8 = sum_with_sign(list8)
    sum9 = sum_with_sign(list9)

    # spectral norm: ||A||_spectral = sqrt[lambda_max(A^dag A)] TODO fix as in 3q_Trotter_compare.py
    print("||sum1||:", np.sqrt(np.max((sum1.dag()*sum1).eigenstates()[0]))) # g12g13
    print("||sum2||:", np.sqrt(np.max((sum2.dag()*sum2).eigenstates()[0]))) # g12g23
    print("||sum3||:", np.sqrt(np.max((sum3.dag()*sum3).eigenstates()[0]))) # g13g23
    print("||sum4||:", np.sqrt(np.max((sum4.dag()*sum4).eigenstates()[0]))) # g12g13
    print("||sum5||:", np.sqrt(np.max((sum5.dag()*sum5).eigenstates()[0]))) # g12g23
    print("||sum6||:", np.sqrt(np.max((sum6.dag()*sum6).eigenstates()[0]))) # g13g23
    print("||sum7||:", np.sqrt(np.max((sum7.dag()*sum7).eigenstates()[0]))) # g12g13
    print("||sum8||:", np.sqrt(np.max((sum8.dag()*sum8).eigenstates()[0]))) # g12g23
    print("||sum9||:", np.sqrt(np.max((sum9.dag()*sum9).eigenstates()[0]))) # g13g23







if __name__ == "__main__":
    main()
