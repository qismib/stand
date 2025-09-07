import numpy as np
import qutip as qt
from itertools import product, permutations

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
    ('X', 'Y'): (+2j, 'Z'),
    ('Y', 'Z'): (+2j, 'X'),
    ('Z', 'X'): (+2j, 'Y'),
    ('Y', 'X'): (-2j, 'Z'),
    ('Z', 'Y'): (-2j, 'X'),
    ('X', 'Z'): (-2j, 'Y'),
    ('0', 'I'): (1, '0'),
    ('0', 'X'): (1, '0'),
    ('0', 'Y'): (1, '0'),
    ('X', '0'): (1, '0'),
    ('I', '0'): (1, '0'),
    ('X', '0'): (1, '0'),
    ('Y', '0'): (1, '0'),
    ('Z', '0'): (1, '0')
}

# ---------------------------------------- PauliCommutators ----------------------------------------

class PauliCommutators():
    dim = 1 # number of qubits in the system, default is 1

    def __init__(self, N):
        self.dim = N

    def generate_pcombo(self, elems=['I','X','Y','Z']):
        """
        `generate_pcombo(self, elems=['I','X','Y','Z'])`
        elems (list): generates all possible combinations of (a subset of) Pauli matrices and the identity.

        product(list, len) generates all possible combinations of length len using elements of list
        """
        return [''.join(p) for p in product(elems, repeat=self.dim)]


    def comm_pstr(self, pstr1, pstr2):
        """
        `comm_pstr(self, pstr1, pstr2)`
        pstr1, pstr2 (tuple): (numeric coefficient, couplings, individual Pauli string)

        comm_pstr function computes the commutator [pstr1, pstr2], keeping track of resulting numeric coefficient and coupling combination.
        TODO: explain how it works.
        return: (numeric coefficient, couplings, individual Pauli string) (tuple)
        """
        coeff1 = pstr1[0]
        coeff2 = pstr2[0]
        ps1 = pstr1[2]
        ps2 = pstr2[2]

        coeff = coeff1 * coeff2
        psfinal = ""
        acomm_cnt = 0

        for i in range(len(ps1)):
            p1 = ps1[i]
            p2 = ps2[i]

            if p1 != p2 and p1 != 'I' and p2 != 'I': acomm_cnt += 1

            pair = pauli_mult_table[(p1, p2)]
            coeff = coeff * pair[0]
            psfinal += pair[1]

        if acomm_cnt % 2 == 0:
            return (1, "", '0'*len(ps1)) 
        
        return (coeff, f"({pstr1[1]})({pstr2[1]})", psfinal)
    
    def comm_lincombo(self, psum1, psum2):
        """
        `comm_lincombo(self, psum1, psum2)`
        psum1, psum2 (list of tuples): [(numeric coefficient 1, couplings 1, individual Pauli string 1), ...]

        comm_lincombo computes the commutator of two sums of Pauli strings, keeping track of numerical coefficients and coupling powers.
        return: [(resulting numeric coefficient 1, resulting couplings 1, resulting individual Pauli string 1), ...] (list of tuples)
        """
        list_comm = []

        for elem1 in psum1:
            for elem2 in psum2:
                res = self.comm_pstr(elem1, elem2)
                list_comm.append((res[0], res[1], res[2]))

        return list_comm

    def refactor_coeff(self, terms):
        """
        `refactor_coeff(self, terms)`
        terms (list of lists): it is a list with elements of the form (numerical coefficient, 'couplings'), which are all the contributions to a certain Pauli string

        return: [(numeric coeff, powers of coupling), ...] list of tuples 
        """
        sum_by_powers = {}

        for coeff, expr in terms:
            idx = expr.replace('(', '').replace(')', '').split("g")[1:]
            powers = (idx.count('12'), idx.count('13'), idx.count('23'))

            if powers in sum_by_powers:
                sum_by_powers[powers] += coeff
            else:
                sum_by_powers[powers] = coeff

        # Convert to desired list format
        summed_terms = [[c, p] for p, c in sum_by_powers.items()]

        def format_term(coeff, powers):
            labels = ['g12', 'g13', 'g23']
            term = ''.join(f'{labels[i]}' * p for i, p in enumerate(powers))
            return f'{coeff}{term if term else ""}'
        
        formatted = [format_term(c, p) for c, p in summed_terms]
        return '(' + ' + '.join(formatted) + ')'


    def nice_print(self, list_terms, title=""):
        """
        `nice_print(self, list_terms, title=None)`
        list_terms (list of tuples): list of tuples (coeff, Pauli string), with coeff and Pauli string both of type str. 
                                     coeff may contain both numerical coefficient and the couplings.
        title (str): string that is printed before list_terms
        """
        
        if title != "":
            print(title)

        for term in list_terms:
                print(f"{self.refactor_coeff(term[0])}{term[1]}")
        
        return
    
    def mulscal_pstr(self, scalar, pstr):
        '''
        implements multiplication by scalar for a pstr (list of tuples (numerical coefficient, couplings, individual Pauli string)) 
        scalar (float)
        '''
        print(pstr[0][0])
        return [(str(float((elem[0][0].real))*scalar), elem[1], elem[2]) for elem in pstr]
    
    def sum_pstr(self, pstr1, pstr2):
        '''
        pstr1, pstr2 (list of tuples): [(numeric coefficient, couplings, individual Pauli string), ...]
        Sums pstr1 and pstr2 taking care of refactoring numerical coefficients.
        return pstr1+pstr2
        '''

        grouped = {}

        for coeff, couplings, pauli in pstr1+pstr2:
            key = (couplings, pauli)
            if key in grouped:
                grouped[key] += coeff
            else:
                grouped[key] = coeff

        return [(coeff, key[0], key[1]) for key, coeff in grouped.items() if coeff != 0]
    


    


# ---------------------------------------- CommutatorNorms ----------------------------------------

class CommutatorNorms():
    N = 1 # number of qubits
    pbasis = {
              'I': qt.qeye(2), 
              'X': qt.sigmax(),
              'Y': qt.sigmay(),
              'Z': qt.sigmaz(),
              '0': qt.Qobj(np.zeros((2, 2), dtype=complex)) # ([0, 0], [0, 0])
             }


    def __init__(self, num_qubits):
        self.N = num_qubits

        
    def pauli_tensor(self, pauli_str):
        """
        Create tensor product of Pauli matrices from string
        s"""
        if len(pauli_str) != self.N:
            raise ValueError(f"String length must match number of qubits ({self.N})")
        
        op = self.pbasis[pauli_str[0]]
        for p in pauli_str[1:]:
            op = qt.tensor(op, self.pbasis[p])
        return op
    
    # ADDITIONAL OBJECTS
    
    def zero_obj(self):
        return self.pauli_tensor('0'*self.N)
    
    def Rx(self, theta):
        return (-1j * theta/2 * self.pbasis['X']).expm()

    def Ry(self, theta):
        return (-1j * theta/2 * self.pbasis['Y']).expm()

    def Rz(self, theta):
        return (-1j * theta/2 * self.pbasis['Z']).expm()
        
    def gen_all_commutators(self, op_list):
        """
        gen_all_commutators(self, op_list)
        Generates all possible commutators given the operators in `op_list`
        """
        dict = {}

        for pair in list(permutations(op_list, 2)):
            first = pair[0]
            secnd = pair[1]
            comm = self.commutator(first, secnd)
            dict[f"[{first},{secnd}]"] = comm

        return dict