from error_value_functions import *

def main1():

    XY3, XZ3, YZ3 = gen_Heisenberg_terms(3)
    XY4, XZ4, YZ4 = gen_Heisenberg_terms(4)
    XY5, XZ5, YZ5 = gen_Heisenberg_terms(5)
    XY6, XZ6, YZ6 = gen_Heisenberg_terms(6)

    XY3XZ3 = comm_lincombo(XY3, XZ3)
    XY3YZ3 = comm_lincombo(XY3, YZ3)
    XZ3YZ3 = comm_lincombo(XZ3, YZ3)

    XY4XZ4 = comm_lincombo(XY4, XZ4)
    XY4YZ4 = comm_lincombo(XY4, YZ4)
    XZ4YZ4 = comm_lincombo(XZ4, YZ4)

    XY5XZ5 = comm_lincombo(XY5, XZ5)
    XY5YZ5 = comm_lincombo(XY5, YZ5)
    XZ5YZ5 = comm_lincombo(XZ5, YZ5)

    XY6XZ6 = comm_lincombo(XY6, XZ6)
    XY6YZ6 = comm_lincombo(XY6, YZ6)
    XZ6YZ6 = comm_lincombo(XZ6, YZ6)


    terms_XY3XZ3 = group_by_pauli_string(XY3XZ3, N=3)
    terms_XY3YZ3 = group_by_pauli_string(XY3YZ3, N=3)
    terms_XZ3YZ3 = group_by_pauli_string(XZ3YZ3, N=3)

    terms_XY4XZ4 = group_by_pauli_string(XY4XZ4, N=4)
    terms_XY4YZ4 = group_by_pauli_string(XY4YZ4, N=4)
    terms_XZ4YZ4 = group_by_pauli_string(XZ4YZ4, N=4)

    terms_XY5XZ5 = group_by_pauli_string(XY5XZ5, N=5)
    terms_XY5YZ5 = group_by_pauli_string(XY5YZ5, N=5)
    terms_XZ5YZ5 = group_by_pauli_string(XZ5YZ5, N=5)

    terms_XY6XZ6 = group_by_pauli_string(XY6XZ6, N=6)
    terms_XY6YZ6 = group_by_pauli_string(XY6YZ6, N=6)
    terms_XZ6YZ6 = group_by_pauli_string(XZ6YZ6, N=6)
    
    mu = 1
    peakness = 0.9
    dt = 0.01


    couplings3, angles3 = fp_couplings(3, peakness, mu)
    couplings4, angles4 = fp_couplings(4, peakness, mu)
    couplings5, angles5 = fp_couplings(5, peakness, mu)
    couplings6, angles6 = fp_couplings(6, peakness, mu)

    # additive error bound

    err_vl_additive = eval_additive_error_term_T1([terms_XY3XZ3, terms_XY3YZ3, terms_XZ3YZ3], couplings3, dt) 
    print(f"additive error value: {err_vl_additive}")

    # tight error bound

    terms_XYXZ_XYYZ3 = merge_terms_dicts([terms_XY3XZ3, terms_XY3YZ3])
    terms_XYXZ_XYYZ3 = merge_terms_dicts([terms_XY3XZ3, terms_XY3YZ3])
    terms_XYXZ_XYYZ3 = merge_terms_dicts([terms_XY3XZ3, terms_XY3YZ3])

    #err_vl_tight = eval_tight_error_term_T1()


def main2():

    XY3, XZ3, YZ3 = gen_Heisenberg_terms(3)
    XY4, XZ4, YZ4 = gen_Heisenberg_terms(4)
    XY5, XZ5, YZ5 = gen_Heisenberg_terms(5)
    XY6, XZ6, YZ6 = gen_Heisenberg_terms(6)

    XY3YZ3 = comm_lincombo(XY3, YZ3)
    YZ3XY3 = comm_lincombo(YZ3, XY3)
    XY3XZ3 = comm_lincombo(XY3, XZ3)
    XZ3XY3 = comm_lincombo(XZ3, XY3)
    XZ3YZ3 = comm_lincombo(XZ3, YZ3)
    YZ3XZ3 = comm_lincombo(YZ3, XZ3)
    

    XY4YZ4 = comm_lincombo(XY4, YZ4)
    YZ4XY4 = comm_lincombo(YZ4, XY4)
    XY4XZ4 = comm_lincombo(XY4, XZ4)
    XZ4XY4 = comm_lincombo(XZ4, XY4)
    XZ4YZ4 = comm_lincombo(XZ4, YZ4)
    YZ4XZ4 = comm_lincombo(YZ4, XZ4)
    
    XY5YZ5 = comm_lincombo(XY5, YZ5)
    YZ5XY5 = comm_lincombo(YZ5, XY5)
    XY5XZ5 = comm_lincombo(XY5, XZ5)
    XZ5XY5 = comm_lincombo(XZ5, XY5)
    XZ5YZ5 = comm_lincombo(XZ5, YZ5)
    YZ5XZ5 = comm_lincombo(YZ5, XZ5)
    
    XY6YZ6 = comm_lincombo(XY6, YZ6)
    YZ6XY6 = comm_lincombo(YZ6, XY6)
    XY6XZ6 = comm_lincombo(XY6, XZ6)
    XZ6XY6 = comm_lincombo(XZ6, XY6)
    XZ6YZ6 = comm_lincombo(XZ6, YZ6)
    YZ6XZ6 = comm_lincombo(YZ6, XZ6)


    terms_XY3YZ3 = group_by_pauli_string(XY3YZ3, N=3)
    terms_YZ3XY3 = group_by_pauli_string(YZ3XY3, N=3)
    terms_XY3XZ3 = group_by_pauli_string(XY3XZ3, N=3)
    terms_XZ3XY3 = group_by_pauli_string(XZ3XY3, N=3)
    terms_XZ3YZ3 = group_by_pauli_string(XZ3YZ3, N=3)
    terms_YZ3XZ3 = group_by_pauli_string(YZ3XZ3, N=3)

    terms_XY4YZ4 = group_by_pauli_string(XY4YZ4, N=4)
    terms_YZ4XY4 = group_by_pauli_string(YZ4XY4, N=4)
    terms_XY4XZ4 = group_by_pauli_string(XY4XZ4, N=4)
    terms_XZ4XY4 = group_by_pauli_string(XZ4XY4, N=4)
    terms_XZ4YZ4 = group_by_pauli_string(XZ4YZ4, N=4)
    terms_YZ4XZ4 = group_by_pauli_string(YZ4XZ4, N=4)

    terms_XY5YZ5 = group_by_pauli_string(XY5YZ5, N=5)
    terms_YZ5XY5 = group_by_pauli_string(YZ5XY5, N=5)
    terms_XY5XZ5 = group_by_pauli_string(XY5XZ5, N=5)
    terms_XZ5XY5 = group_by_pauli_string(XZ5XY5, N=5)
    terms_XZ5YZ5 = group_by_pauli_string(XZ5YZ5, N=5)
    terms_YZ5XZ5 = group_by_pauli_string(YZ5XZ5, N=5)

    terms_XY6YZ6 = group_by_pauli_string(XY6YZ6, N=6)
    terms_YZ6XY6 = group_by_pauli_string(YZ6XY6, N=6)
    terms_XY6XZ6 = group_by_pauli_string(XY6XZ6, N=6)
    terms_XZ6XY6 = group_by_pauli_string(XZ6XY6, N=6)
    terms_XZ6YZ6 = group_by_pauli_string(XZ6YZ6, N=6)
    terms_YZ6XZ6 = group_by_pauli_string(YZ6XZ6, N=6)

    # tight error bound

    terms_XZXY_YZXY3 = merge_terms_dicts([terms_XZ3XY3, terms_YZ3XY3])
    terms_XYXZ_YZXZ3 = merge_terms_dicts([terms_XY3XZ3, terms_YZ3XZ3])
    terms_XYYZ_XZYZ3 = merge_terms_dicts([terms_XY3YZ3, terms_XZ3YZ3])

    dt = 0.01
    mu = 1
    peakness = 0.9

    couplings3, angles3 = fp_couplings(3, peakness, mu)

    err_vl = eval_tight_error_term_T1([terms_XZXY_YZXY3, terms_XYXZ_YZXZ3, terms_XYYZ_XZYZ3], couplings3, dt)
    print(f"tight error value: {err_vl}")

if __name__ == "__main__":
    main1()
    main2()














