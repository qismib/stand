# ---------------------------- T2_commutators_to_latex.py ---------------------------
#------------------------------------------------------------------------------------
# This scripts computes second order Trotter commutators and formats the result into 
# a (almost) LaTeX compatible string.

# I aim to do N = 3,4,5,6 simulations for the first order Trotter error.
# Note: N=2 additive Trotter error is zero.

import numpy as np
import qutip as qt
import math
import re
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
        s = str(s).strip("(").strip(")")
        is_complex = False
        re_part = "0"
        im_part = "0"

        def simplify_decimal(d):
            int_d, dec_d = d.split(".")
            num_d = int(int_d + dec_d)
            den_d = 10 ** len(dec_d) 
            g_d = math.gcd(num_d, den_d)
            num_d //= g_d
            den_d //= g_d 
            return str(num_d), str(den_d)

        if 'j' in s: 
            match_im = re.search(r'([-+]?\d*\.?\d+)(?=j)', s)
            match_re = re.search(r'^([-+]?\d*\.?\d+)(?=[+-]\d*\.?\d*j$)', s)
            if match_im and float(match_im.group(1)) != 0 and match_re:
                re_part = match_re.group(1)
                im_part = match_im.group(1)
                is_complex = True
            elif match_im and float(match_im.group(1)) != 0:
                re_part = "0"
                im_part = match_im.group(1)
                is_complex = True
            elif match_re:
                re_part = match_re.group(1)
                im_part = "0"
        else:
            im_part = 0
            re_part = s
            if "." in re_part: 
                num_re, den_re = simplify_decimal(re_part)
                return f"(\\frac{{{num_re}}}{{{den_re}}} + 0i)"
            else:
                return f"({re_part} + 0i)"

        if re_part == "0" and im_part == "0": return f"0"

        if "." not in re_part and "." not in im_part:
            if re_part == "0" and not is_complex:
                return "0"
            elif not is_complex:
                return f"({re_part}+0i)"
            else:
                return f"({re_part}+{im_part}i)"
        elif "." not in re_part:
            if is_complex:
                num_im, den_im = simplify_decimal(im_part)
                if "-" in num_im:
                    return f"({re_part}+\\frac{{{num_im}}}{{{den_im}}}i)"
                else:
                    return f"({re_part}+\\frac{{+{num_im}}}{{{den_im}}}i)"
            else:
                return f"({re_part}+0i)"
        elif "." not in im_part:
            if is_complex:
                num_re, den_re = simplify_decimal(re_part)
                if "-" in num_re:
                    return f"(\\frac{{{num_re}}}{{{den_re}}}+{im_part}i)"
                else:
                    return f"(\\frac{{+{num_re}}}{{{den_re}}}+{im_part}i)"
            else:
                return f"({re_part}+0j)"
        else:
            if is_complex:
                num_re, den_re = simplify_decimal(re_part)
                num_im, den_im = simplify_decimal(im_part)
                if "-" in num_re and "-" in num_im:
                    return f"(\\frac{{{num_re}}}{{{den_re}}}+\\frac{{{num_im}}}{{{den_im}}}i)"
                elif "-" in num_re:
                    return f"(\\frac{{{num_re}}}{{{den_re}}}+\\frac{{+{num_im}}}{{{den_im}}}i)"
                elif "-" in num_im:
                    return f"\\frac{{+{num_re}}}{{{den_re}}}+\\frac{{{num_im}}}{{{den_im}}}i"
                else:
                    return f"(\\frac{{+{num_re}}}{{{den_re}}}+\\frac{{+{num_im}}}{{{den_im}}}i)"
            else:
                if "-" in num_re:
                    return f"(\\frac{{{num_re}}}{{{den_re}}}+0i)"
                else:
                    return f"(\\frac{{+{num_re}}}{{{den_re}}}+0i)"

        # NOTE: I deal with signs in this (apparently useless) way in order to get a consistent output
        # (can use multicursor easily)
            
    # -------------------------------------------------------------------------------------------
    # actual function
    latex_str = f"{title}="
    for key, value in terms.items():
        temp = "\\left|"
        for elem in value:
            ref_coeff = clean_coeff(elem[1])
            temp += str_to_latex_frac(elem[0])
            if "(" not in str_to_latex_frac(elem[0]) or ")" not in str_to_latex_frac(elem[0]): raise Exception(f"something wrong: {str_to_latex_frac(elem[0])}")
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


    terms_B3_C3A3 = group_by_pauli_string(B3_C3A3, N=3)
    terms_B3_B3A3 = group_by_pauli_string(B3_B3A3, N=3)
    terms_C3_C3A3 = group_by_pauli_string(C3_C3A3, N=3)
    terms_C3_B3A3 = group_by_pauli_string(C3_B3A3, N=3)
    terms_A3_A3B3 = group_by_pauli_string(A3_A3B3, N=3)
    terms_C3_C3B3 = group_by_pauli_string(C3_C3B3, N=3)
    terms_B3_B3C3 = group_by_pauli_string(B3_B3C3, N=3)
    terms_A3_A3C3 = group_by_pauli_string(A3_A3C3, N=3)

    terms_B4_C4A4 = group_by_pauli_string(B4_C4A4, N=4)
    terms_B4_B4A4 = group_by_pauli_string(B4_B4A4, N=4)
    terms_C4_C4A4 = group_by_pauli_string(C4_C4A4, N=4)
    terms_C4_B4A4 = group_by_pauli_string(C4_B4A4, N=4)
    terms_A4_A4B4 = group_by_pauli_string(A4_A4B4, N=4)
    terms_C4_C4B4 = group_by_pauli_string(C4_C4B4, N=4)
    terms_B4_B4C4 = group_by_pauli_string(B4_B4C4, N=4)
    terms_A4_A4C4 = group_by_pauli_string(A4_A4C4, N=4)

    terms_B5_C5A5 = group_by_pauli_string(B5_C5A5, N=5)
    terms_B5_B5A5 = group_by_pauli_string(B5_B5A5, N=5)
    terms_C5_C5A5 = group_by_pauli_string(C5_C5A5, N=5)
    terms_C5_B5A5 = group_by_pauli_string(C5_B5A5, N=5)
    terms_A5_A5B5 = group_by_pauli_string(A5_A5B5, N=5)
    terms_C5_C5B5 = group_by_pauli_string(C5_C5B5, N=5)
    terms_B5_B5C5 = group_by_pauli_string(B5_B5C5, N=5)
    terms_A5_A5C5 = group_by_pauli_string(A5_A5C5, N=5)

    terms_B6_C6A6 = group_by_pauli_string(B6_C6A6, N=6)
    terms_B6_B6A6 = group_by_pauli_string(B6_B6A6, N=6)
    terms_C6_C6A6 = group_by_pauli_string(C6_C6A6, N=6)
    terms_C6_B6A6 = group_by_pauli_string(C6_B6A6, N=6)
    terms_A6_A6B6 = group_by_pauli_string(A6_A6B6, N=6)
    terms_C6_C6B6 = group_by_pauli_string(C6_C6B6, N=6)
    terms_B6_B6C6 = group_by_pauli_string(B6_B6C6, N=6)
    terms_A6_A6C6 = group_by_pauli_string(A6_A6C6, N=6)

    latex_B3_C3A3 = dict_to_latex(terms_B3_C3A3, "[B_3,[C_3,A_3]]")
    latex_B3_B3A3 = dict_to_latex(terms_B3_B3A3, "[B_3,[B_3,A_3]]")
    latex_C3_C3A3 = dict_to_latex(terms_C3_C3A3, "[C_3,[C_3,A_3]]")
    latex_C3_B3A3 = dict_to_latex(terms_C3_B3A3, "[C_3,[B_3,A_3]]")
    latex_A3_A3B3 = dict_to_latex(terms_A3_A3B3, "[A_3,[A_3,B_3]]")
    latex_C3_C3B3 = dict_to_latex(terms_C3_C3B3, "[C_3,[C_3,B_3]]")
    latex_B3_B3C3 = dict_to_latex(terms_B3_B3C3, "[B_3,[B_3,C_3]]")
    latex_A3_A3C3 = dict_to_latex(terms_A3_A3C3, "[A_3,[A_3,C_3]]")
    latex_B4_C4A4 = dict_to_latex(terms_B4_C4A4, "[B_4,[C_4,A_4]]")
    latex_B4_B4A4 = dict_to_latex(terms_B4_B4A4, "[B_4,[B_4,A_4]]")
    latex_C4_C4A4 = dict_to_latex(terms_C4_C4A4, "[C_4,[C_4,A_4]]")
    latex_C4_B4A4 = dict_to_latex(terms_C4_B4A4, "[C_4,[B_4,A_4]]")
    latex_A4_A4B4 = dict_to_latex(terms_A4_A4B4, "[A_4,[A_4,B_4]]")
    latex_C4_C4B4 = dict_to_latex(terms_C4_C4B4, "[C_4,[C_4,B_4]]")
    latex_B4_B4C4 = dict_to_latex(terms_B4_B4C4, "[B_4,[B_4,C_4]]")
    latex_A4_A4C4 = dict_to_latex(terms_A4_A4C4, "[A_4,[A_4,C_4]]")
    latex_B5_C5A5 = dict_to_latex(terms_B5_C5A5, "[B_5,[C_5,A_5]]")
    latex_B5_B5A5 = dict_to_latex(terms_B5_B5A5, "[B_5,[B_5,A_5]]")
    latex_C5_C5A5 = dict_to_latex(terms_C5_C5A5, "[C_5,[C_5,A_5]]")
    latex_C5_B5A5 = dict_to_latex(terms_C5_B5A5, "[C_5,[B_5,A_5]]")
    latex_A5_A5B5 = dict_to_latex(terms_A5_A5B5, "[A_5,[A_5,B_5]]")
    latex_C5_C5B5 = dict_to_latex(terms_C5_C5B5, "[C_5,[C_5,B_5]]")
    latex_B5_B5C5 = dict_to_latex(terms_B5_B5C5, "[B_5,[B_5,C_5]]")
    latex_A5_A5C5 = dict_to_latex(terms_A5_A5C5, "[A_5,[A_5,C_5]]")
    latex_B6_C6A6 = dict_to_latex(terms_B6_C6A6, "[B_6,[C_6,A_6]]")
    latex_B6_B6A6 = dict_to_latex(terms_B6_B6A6, "[B_6,[B_6,A_6]]")
    latex_C6_C6A6 = dict_to_latex(terms_C6_C6A6, "[C_6,[C_6,A_6]]")
    latex_C6_B6A6 = dict_to_latex(terms_C6_B6A6, "[C_6,[B_6,A_6]]")
    latex_A6_A6B6 = dict_to_latex(terms_A6_A6B6, "[A_6,[A_6,B_6]]")
    latex_C6_C6B6 = dict_to_latex(terms_C6_C6B6, "[C_6,[C_6,B_6]]")
    latex_B6_B6C6 = dict_to_latex(terms_B6_B6C6, "[B_6,[B_6,C_6]]")
    latex_A6_A6C6 = dict_to_latex(terms_A6_A6C6, "[A_6,[A_6,C_6]]")
    
        
    print(latex_B4_B4A4)
    


if __name__ == "__main__":
    main()