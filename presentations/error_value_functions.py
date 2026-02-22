import numpy as np
import warnings

def gen_Heisenberg_terms(N):
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


pauli_mult_table = {
    ('I', 'I'): (1, 'I'),
    ('I', 'X'): (1, 'X'),
    ('I', 'Y'): (1, 'Y'),
    ('I', 'Z'): (1, 'Z'),
    ('X', 'I'): (1, 'X'),
    ('Y', 'I'): (1, 'Y'),
    ('Z', 'I'): (1, 'Z'),
    ('X', 'X'): (1, 'I'),
    ('Y', 'Y'): (1, 'I'),
    ('Z', 'Z'): (1, 'I'),
    ('Y', 'Z'): (+1j, 'X'),
    ('X', 'Y'): (+1j, 'Z'),
    ('Z', 'X'): (+1j, 'Y'),
    ('Y', 'X'): (-1j, 'Z'),
    ('Z', 'Y'): (-1j, 'X'),
    ('X', 'Z'): (-1j, 'Y'),
    ('0', 'I'): (0, '0'),
    ('0', 'X'): (0, '0'),
    ('0', 'Y'): (0, '0'),
    ('X', '0'): (0, '0'),
    ('I', '0'): (0, '0'),
    ('X', '0'): (0, '0'),
    ('Y', '0'): (0, '0'),
    ('Z', '0'): (0, '0')
}


def comm_pstr(pstr1, pstr2):
    """
    `comm_pstr(self, pstr1, pstr2)`
    pstr1, pstr2 (tuple): (numeric coefficient, couplings, individual Pauli string)

    comm_pstr function computes the commutator [pstr1, pstr2], keeping track of resulting numeric coefficient and coupling combination.
    return: (numeric coefficient, couplings, individual Pauli string) (tuple)
    """
    ps1 = pstr1[2] # e.g. "IXX"
    ps2 = pstr2[2] # e.g. "YIY"

    coeff = pstr1[0] * pstr2[0]
    psfinal = ""
    acomm_cnt = 0

    for i in range(len(ps1)):
        p1 = ps1[i]
        p2 = ps2[i]

        if p1 != p2 and p1 != 'I' and p2 != 'I': acomm_cnt += 1

        pair = pauli_mult_table[(p1, p2)]
        coeff = coeff * pair[0]
        psfinal += pair[1]

    if acomm_cnt % 2 == 0 or '0' in psfinal:
        return (0, "", '0'*len(ps1))
    
    return (2*coeff, f"({pstr1[1]})({pstr2[1]})", psfinal) # factor of 2 is introduced here


def comm_lincombo(psum1, psum2):
    """
    `comm_lincombo(self, psum1, psum2)`
    psum1, psum2 (list of tuples): [(numeric coefficient 1, couplings 1, individual Pauli string 1), ...]

    comm_lincombo computes the commutator of two sums of Pauli strings, keeping track of numerical coefficients and coupling powers.
    return: [(resulting numeric coefficient 1, resulting couplings 1, resulting individual Pauli string 1), ...] (list of tuples)
    """
    list_comm = []

    for elem1 in psum1:
        for elem2 in psum2:
            res = comm_pstr(elem1, elem2)
            if '0' not in res[2]: list_comm.append((res[0], res[1], res[2]))

    return list_comm


def group_by_pauli_string(comm, N):
    """
    `group_by_pauli_string(comm, N)`
    comm (list of tuples): [(numeric coefficient, couplings, Pauli string), ...]
    N (int)              : number of qubits

    group_by_pauli_string groups terms with the same Pauli string, and within those terms sums
    contributions with the same coupling coefficient

    return: dict {Pauli string: [[coefficient, coupling], ...], ...} 
    """

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
                if cnt_not == len(list_terms[term[2]]): list_terms[term[2]].append([term[0], term[1]]) # TODO can it be changed to ([term[0], term[1])?

    return list_terms

PAULI = {
    'I': np.eye(2, dtype=complex),
    'X': np.array([[0, 1], [1, 0]], dtype=complex),
    'Y': np.array([[0, -1j], [1j, 0]], dtype=complex),
    'Z': np.array([[1, 0], [0, -1]], dtype=complex),
}

def pauli_string_to_matrix(pauli_string):
    """
    Convert a Pauli string (e.g. 'XYZ') into its matrix representation.
    Invalid input returns a zero matrix and raises a warning.
    """
    if pauli_string == '':
        warnings.warn("Empty Pauli string provided.")
        return np.zeros((2 ** len(pauli_string), 2 ** len(pauli_string)), dtype=complex)
    
    for p in pauli_string:
        if p not in PAULI:
            warnings.warn(f"Invalid Pauli operator '{p}' in string '{pauli_string}'.")
            return np.zeros((2 ** len(pauli_string), 2 ** len(pauli_string)), dtype=complex)

    matrix = PAULI[pauli_string[0]]
    for p in pauli_string[1:]:
        matrix = np.kron(matrix, PAULI[p])

    return matrix

def eval_additive_error_term_T1(error_terms, couplings, t): 
    """
    ### Formula for a_1(N,G)t^2

    error_terms (list of dicts): [{"XYZ": [[-0.5j, "(g12)(g23)"], ...], ...}, {"XYZ": [[0.5j, "(g12)(g13)"], ...], ...}]
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

    
    err_vl = 0

    for spec_norm_term in error_terms:
        error_contribution = 0
        for pstring, terms in spec_norm_term.items():
            coeff = 0
            for term in terms:
                # if term[1] != '': coeff += float(str(term[0]).strip('j')) * eval_couplings(term[1], couplings)
                if term[1] != '': coeff += term[0] * eval_couplings(term[1], couplings)
            error_contribution += coeff*pauli_string_to_matrix(pstring)
        err_vl += np.sqrt(np.max(np.linalg.svd(error_contribution, compute_uv=False)))

    err_vl = 0.5 * err_vl * t ** 2

    return err_vl

def merge_terms_dicts(dict_list):
    def canonicalize(term):
        factors = term.replace(")(", " ").replace("(", "").replace(")", "").split()
        factors = sorted(factors)
        return "(" + ")(".join(factors) + ")"

    result = {}

    for d in dict_list:
        for key, pairs in d.items():
            if key not in result:
                result[key] = {}

            for coeff, term in pairs:
                norm_term = canonicalize(term)
                if norm_term in result[key]:
                    result[key][norm_term] += coeff
                else:
                    result[key][norm_term] = coeff

    for key in result:
        result[key] = [[coeff, term] for term, coeff in result[key].items() if coeff != 0]

    return result

def eval_tight_error_term_T1(error_terms, couplings, t): 
    """
    ### Formula for e_1(N,G,t)

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
    
    err_vl = 0

    for spec_norm_term in error_terms:
        error_contribution = 0
        for pstring, terms in spec_norm_term.items():
            coeff = 0
            for term in terms:
                # if term[1] != '': coeff += float(str(term[0]).strip('j')) * eval_couplings(term[1], couplings)
                if term[1] != '': coeff += term[0] * eval_couplings(term[1], couplings)

            error_contribution += coeff*pauli_string_to_matrix(pstring)
        err_vl += np.max(np.linalg.svd(error_contribution, compute_uv=False))

    err_vl = err_vl * t ** 2 / 8

    return err_vl

def theta_nu(N, i, j, vl):
    return np.arccos(vl) * abs(i-j) / (N-1)

def fp_couplings(N, peakness, mu): # fp = forward peaked
    """
    N (int): number of neutrinos
    peakness (float): in [-1,1)
    """
    if peakness < -1 or peakness >= 1: raise Exception(f"Value of peakness argument must be in [-1,1) interval, but {peakness} ins't.")
    coupl_dic = {}
    angle_dic = {}
    for i in range(N):
        for j in range(i+1, N):
            angle_dic[f"g{i+1}{j+1}"] = theta_nu(N, i, j, peakness)
            coupl_dic[f"g{i+1}{j+1}"] = mu * (1-np.cos(theta_nu(N, i, j, peakness))) / N
    
    return coupl_dic, angle_dic

def Theta(angles):
    return max(1 - np.cos(angle) for angle in angles.values())








