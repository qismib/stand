import sys; sys.path.append("../classes")
from CommutatorTool import PauliCommutators


def main():

    ctl3 = PauliCommutators(3)

    A3 = [(0.5, "g12","XXI"), (0.5, "g12", "YYI"), (0.5, "g13", "XIX"), (0.5, "g13", "YIY"), (0.5, "g23", "IXX"), (0.5, "g23", "IYY")]
    B3 = [(0.5, "g12","XXI"), (0.5, "g12", "ZZI"), (0.5, "g13", "XIX"), (0.5, "g13", "ZIZ"), (0.5, "g23", "IXX"), (0.5, "g23", "IZZ")]
    C3 = [(0.5, "g12","YYI"), (0.5, "g12", "ZZI"), (0.5, "g13", "YIY"), (0.5, "g13", "ZIZ"), (0.5, "g23", "IYY"), (0.5, "g23", "IZZ")]

    # first order commutators
    A3B3 = ctl3.comm_lincombo(A3, B3)
    A3C3 = ctl3.comm_lincombo(A3, C3)
    B3C3 = ctl3.comm_lincombo(B3, C3)

    # second order commutators
    A3_A3B3 = ctl3.comm_lincombo(A3, A3B3)
    B3_A3B3 = ctl3.comm_lincombo(B3, A3B3)
    C3_A3B3 = ctl3.comm_lincombo(C3, A3B3)
    A3_A3C3 = ctl3.comm_lincombo(A3, A3C3)
    B3_A3C3 = ctl3.comm_lincombo(B3, A3C3)
    C3_A3C3 = ctl3.comm_lincombo(C3, A3C3)
    A3_B3C3 = ctl3.comm_lincombo(A3, B3C3)
    B3_B3C3 = ctl3.comm_lincombo(B3, B3C3)
    C3_B3C3 = ctl3.comm_lincombo(C3, B3C3)

    coeff_XXI = []
    coeff_XIX = []
    coeff_IXX = []
    coeff_YYI = []
    coeff_YIY = []
    coeff_IYY = []
    coeff_ZZI = []
    coeff_ZIZ = []
    coeff_IZZ = []

    zero_cnt = 0

    for elem in A3_A3B3+B3_A3B3+C3_A3B3+A3_A3C3+B3_A3C3+C3_A3C3+A3_B3C3+B3_B3C3+C3_B3C3:
        if '0' not in elem[2]:
            if elem[2] == "XXI": 
                coeff_XXI.append((elem[0], elem[1]))
            elif elem[2] == "XIX": 
                coeff_XIX.append((elem[0], elem[1]))
            elif elem[2] == "IXX": 
                coeff_IXX.append((elem[0], elem[1]))
            elif elem[2] == "YYI": 
                coeff_YYI.append((elem[0], elem[1]))
            elif elem[2] == "YIY": 
                coeff_YIY.append((elem[0], elem[1]))
            elif elem[2] == "IYY": 
                coeff_IYY.append((elem[0], elem[1]))
            elif elem[2] == "ZZI": 
                coeff_ZZI.append((elem[0], elem[1]))
            elif elem[2] == "ZIZ": 
                coeff_ZIZ.append((elem[0], elem[1]))
            elif elem[2] == "IZZ": 
                coeff_IZZ.append((elem[0], elem[1]))
            else: raise Exception(f"ERROR: {elem[2]} is not one of the expected contributions.")
        else:
            zero_cnt += 1

    coeffs = []
    coeffs.append((coeff_XXI, "XXI"))
    coeffs.append((coeff_XIX, "XIX"))
    coeffs.append((coeff_IXX, "IXX"))
    coeffs.append((coeff_YYI, "YYI"))
    coeffs.append((coeff_YIY, "YIY"))
    coeffs.append((coeff_IYY, "IYY"))
    coeffs.append((coeff_ZZI, "ZZI"))
    coeffs.append((coeff_ZIZ, "ZIZ"))
    coeffs.append((coeff_IZZ, "IZZ"))

    
    ctl3.nice_print(coeffs, title="[A3,[A3,B3]]+[B3,[A3,B3]]+[C3,[A3,B3]]+[A3,[A3,C3]]+[B3,[A3,C3]]+[C3,[A3,C3]]+[A3,[B3,C3]]+[B3,[B3,C3]]+[C3,[B3,C3]]:")



if __name__ == "__main__":
    main()
