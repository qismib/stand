# ---------------------------- T1_commutators_to_latex.py ------------------------------
#---------------------------------------------------------------------------------------
# Providing bounds for the additive Trotter error requires computing spectral norms
# of commutators of complicated Hamiltonian operators. 
# In the digital-analog setup I consider, the Hamiltonian is made up of 
# H_{XY}^N, H_{XZ}^N, H_{YZ}^N operators (see Section 5). 
# I'm interested in providing bounds for the scaling of the additive Trotter
# error as a function of the qubit number N. 
# This scripts computes the commutators and formats the result into a (almost) LaTeX
# compatible string.

# I aim to do N = 3,4,5,6 simulations for the first order Trotter error.
# Note: N=2 additive Trotter error is zero.

import numpy as np
import qutip as qt
import math
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg')
import sys; sys.path.append("../classes")
from CommutatorTool import PauliCommutators

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
            for j in range(i + 1, N):
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

# ---------- dict_to_latex(terms, title) ----------

def dict_to_latex(terms, title):
    """
    terms(dict): {"XYZ": [[-0.5j, "(g12)(g23)"], ...], ...}
    """

    # auxiliary functions
    def clean_coeff(coeff_str):
        split_str = coeff_str.split(")(")
        split_str = [p.strip(')').strip('(').strip('g') for p in split_str]

        terms = []
        for term in set(split_str):
            count = sum(1 for pow in split_str if pow == term)
            terms.append((term, count))
            
        return terms
    
    def pstr_to_subscript(pstr):
        substr = "\\right|"
        for i, op in enumerate(pstr):
            if op != 'I': substr += f"{op}_{i+1}"

        return f"{substr}+"
    
    def str_to_latex_frac(s):
        is_complex = False
        s = str(s)
        print(s)

        if 'j' in s: 
            s = s.strip('j')
            is_complex = True

        if "." not in s: return f"{s}/1"

        int_part, dec_part = s.split(".")
        numerator = int(int_part + dec_part)
        denominator = 10 ** len(dec_part)

        g = math.gcd(numerator, denominator)
        numerator //= g
        denominator //= g  

        if is_complex:
            temp = f"\\frac{{{numerator}}}{{{denominator}}}i"
            if "-" in temp:
                return temp
            else:
                return f"\\frac{{+{numerator}}}{{{denominator}}}i"
        else:
            temp = f"\\frac{{{numerator}}}{{{denominator}}}"
            if "-" in temp:
                return temp
            else:
                return f"\\frac{{+{numerator}}}{{{denominator}}}"
            
    # -------------------------------------------------------------------------------------------
    # actual function
    latex_str = f"{title}="
    for key, value in terms.items():
        temp = "\\left|"
        for elem in value:
            ref_coeff = clean_coeff(elem[1])
            temp += str_to_latex_frac(elem[0])
            for coeff in ref_coeff:
                if coeff[1] == 1:
                    temp += f"g_{{{coeff[0]}}}"
                else:
                    temp += f"g_{{{coeff[0]}}}^{coeff[1]}"

            if elem is not value[-1]: temp += "+" 


        temp += pstr_to_subscript(key) 
        latex_str += temp + "\n"

    return latex_str

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


    



def main():

    ctl3 = PauliCommutators(N=3)
    ctl4 = PauliCommutators(N=4)
    ctl5 = PauliCommutators(N=5)
    ctl6 = PauliCommutators(N=6)

    A3, B3, C3 = gen_Heisenberg_terms(3)
    A4, B4, C4 = gen_Heisenberg_terms(4)
    A5, B5, C5 = gen_Heisenberg_terms(5)
    A6, B6, C6 = gen_Heisenberg_terms(6)

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


    latex_A3B3 = dict_to_latex(terms_A3B3, "[A_3,B_3]")
    latex_A3C3 = dict_to_latex(terms_A3C3, "[A_3,C_3]")
    latex_B3C3 = dict_to_latex(terms_B3C3, "[B_3,C_3]")
    latex_A4B4 = dict_to_latex(terms_A4B4, "[A_4,B_4]")
    latex_A4C4 = dict_to_latex(terms_A4C4, "[A_4,C_4]")
    latex_B4C4 = dict_to_latex(terms_B4C4, "[B_4,C_4]")
    latex_A5B5 = dict_to_latex(terms_A5B5, "[A_5,B_5]")
    latex_A5C5 = dict_to_latex(terms_A5C5, "[A_5,C_5]")
    latex_B5C5 = dict_to_latex(terms_B5C5, "[B_5,C_5]")
    latex_A6B6 = dict_to_latex(terms_A6B6, "[A_6,B_6]")
    latex_A6C6 = dict_to_latex(terms_A6C6, "[A_6,C_6]")
    latex_B6C6 = dict_to_latex(terms_B6C6, "[B_6,C_6]")
        
    print(latex_A3B3)
    


if __name__ == "__main__":
    main()