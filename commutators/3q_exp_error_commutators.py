# --------------------------- 3q_exp_error_commutators.py --------------------------
# ----------------------------------------------------------------------------------
# Compute the commutators that appear in section 4.2.2 relevant to the second order
# exponentiated error estimate.

import sys; sys.path.append("../classes")
from CommutatorTool import PauliCommutators


def main():

    ctl3 = PauliCommutators(3)

    A3 = [(0.5, "g12","XXI"), (0.5, "g12", "YYI"), (0.5, "g13", "XIX"), (0.5, "g13", "YIY"), (0.5, "g23", "IXX"), (0.5, "g23", "IYY")]
    B3 = [(0.5, "g12","XXI"), (0.5, "g12", "ZZI"), (0.5, "g13", "XIX"), (0.5, "g13", "ZIZ"), (0.5, "g23", "IXX"), (0.5, "g23", "IZZ")]
    C3 = [(0.5, "g12","YYI"), (0.5, "g12", "ZZI"), (0.5, "g13", "YIY"), (0.5, "g13", "ZIZ"), (0.5, "g23", "IYY"), (0.5, "g23", "IZZ")]

    # first order commutators
    B3A3 = ctl3.comm_lincombo(B3, A3)
    C3A3 = ctl3.comm_lincombo(C3, A3)
    C3B3 = ctl3.comm_lincombo(C3, B3)
    B3C3 = ctl3.comm_lincombo(B3, C3)
    A3B3 = ctl3.comm_lincombo(A3, B3)
    A3C3 = ctl3.comm_lincombo(A3, C3)

    # second order commutators

    B3_B3A3 = ctl3.comm_lincombo(B3, B3A3)
    B3_C3A3 = ctl3.comm_lincombo(B3, C3A3)
    C3_B3A3 = ctl3.comm_lincombo(C3, B3A3)
    C3_C3A3 = ctl3.comm_lincombo(C3, C3A3)
    C3_C3B3 = ctl3.comm_lincombo(C3, C3B3)
    B3_B3C3 = ctl3.comm_lincombo(B3, B3C3)
    A3_A3B3 = ctl3.comm_lincombo(A3, A3B3)
    A3_A3C3 = ctl3.comm_lincombo(A3, A3C3)
    

    # Expression I want to compute: 1/12 ([C,[C,B]]+[B,[B,A]]+[B,[C,A]]+[C,[B,A]]+[C,[C,A]]) - 1/24 ([B,[B,C]]+[A,[A,B]]+[A,[A,C]]) = U/12 - D/24

    # compute U 
    coeffU_XXI = []
    coeffU_XIX = []
    coeffU_IXX = []
    coeffU_YYI = []
    coeffU_YIY = []
    coeffU_IYY = []
    coeffU_ZZI = []
    coeffU_ZIZ = []
    coeffU_IZZ = []

    zero_cnt = 0

    for elem in C3_C3B3+B3_B3A3+B3_C3A3+C3_B3A3+C3_C3A3:
        if '0' not in elem[2]:
            if elem[2] == "XXI": 
                coeffU_XXI.append((elem[0], elem[1]))
            elif elem[2] == "XIX": 
                coeffU_XIX.append((elem[0], elem[1]))
            elif elem[2] == "IXX": 
                coeffU_IXX.append((elem[0], elem[1]))
            elif elem[2] == "YYI": 
                coeffU_YYI.append((elem[0], elem[1]))
            elif elem[2] == "YIY": 
                coeffU_YIY.append((elem[0], elem[1]))
            elif elem[2] == "IYY": 
                coeffU_IYY.append((elem[0], elem[1]))
            elif elem[2] == "ZZI": 
                coeffU_ZZI.append((elem[0], elem[1]))
            elif elem[2] == "ZIZ": 
                coeffU_ZIZ.append((elem[0], elem[1]))
            elif elem[2] == "IZZ": 
                coeffU_IZZ.append((elem[0], elem[1]))
            else: raise Exception(f"ERROR: {elem[2]} is not one of the expected contributions.")
        else:
            zero_cnt += 1

    coeffsU = []
    coeffsU.append((coeffU_XXI, "XXI"))
    coeffsU.append((coeffU_XIX, "XIX"))
    coeffsU.append((coeffU_IXX, "IXX"))
    coeffsU.append((coeffU_YYI, "YYI"))
    coeffsU.append((coeffU_YIY, "YIY"))
    coeffsU.append((coeffU_IYY, "IYY"))
    coeffsU.append((coeffU_ZZI, "ZZI"))
    coeffsU.append((coeffU_ZIZ, "ZIZ"))
    coeffsU.append((coeffU_IZZ, "IZZ"))
    
    #ctl3.nice_print(coeffsU, title="[A3,[A3,B3]]+[B3,[A3,B3]]+[C3,[A3,B3]]+[A3,[A3,C3]]+[B3,[A3,C3]]+[C3,[A3,C3]]+[A3,[B3,C3]]+[B3,[B3,C3]]+[C3,[B3,C3]]:")

    # compute D 
    coeffV_XXI = []
    coeffV_XIX = []
    coeffV_IXX = []
    coeffV_YYI = []
    coeffV_YIY = []
    coeffV_IYY = []
    coeffV_ZZI = []
    coeffV_ZIZ = []
    coeffV_IZZ = []

    zero_cnt = 0

    for elem in B3_B3C3+A3_A3B3+A3_A3C3:
        if '0' not in elem[2]:
            if elem[2] == "XXI": 
                coeffV_XXI.append((elem[0], elem[1]))
            elif elem[2] == "XIX": 
                coeffV_XIX.append((elem[0], elem[1]))
            elif elem[2] == "IXX": 
                coeffV_IXX.append((elem[0], elem[1]))
            elif elem[2] == "YYI": 
                coeffV_YYI.append((elem[0], elem[1]))
            elif elem[2] == "YIY": 
                coeffV_YIY.append((elem[0], elem[1]))
            elif elem[2] == "IYY": 
                coeffV_IYY.append((elem[0], elem[1]))
            elif elem[2] == "ZZI": 
                coeffV_ZZI.append((elem[0], elem[1]))
            elif elem[2] == "ZIZ": 
                coeffV_ZIZ.append((elem[0], elem[1]))
            elif elem[2] == "IZZ": 
                coeffV_IZZ.append((elem[0], elem[1]))
            else: raise Exception(f"ERROR: {elem[2]} is not one of the expected contributions.")
        else:
            zero_cnt += 1

    coeffsV = []
    coeffsV.append((coeffV_XXI, "XXI"))
    coeffsV.append((coeffV_XIX, "XIX"))
    coeffsV.append((coeffV_IXX, "IXX"))
    coeffsV.append((coeffV_YYI, "YYI"))
    coeffsV.append((coeffV_YIY, "YIY"))
    coeffsV.append((coeffV_IYY, "IYY"))
    coeffsV.append((coeffV_ZZI, "ZZI"))
    coeffsV.append((coeffV_ZIZ, "ZIZ"))
    coeffsV.append((coeffV_IZZ, "IZZ"))

    print("Second order contributions:")
    ctl3.nice_print(coeffsU, "U term: (coeff 1/12)")
    ctl3.nice_print(coeffsV, "V term: (coeff -1/24)")

    # TODO: - implement scalar multiplication
    #       - compute directly 1/12 U - 1/24 V with all the necessary refactoring
    



if __name__ == "__main__":
    main()
