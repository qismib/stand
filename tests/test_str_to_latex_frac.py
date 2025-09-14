
import re
import math

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
            return f"{re_part}"
        else:
            return f"{re_part}+{im_part}i"
    elif "." not in re_part:
        if is_complex:
            num_im, den_im = simplify_decimal(im_part)
            if "-" in num_im:
                return f"({re_part}+\\frac{{{num_im}}}{{{den_im}}}i)"
            else:
                return f"({re_part}+\\frac{{+{num_im}}}{{{den_im}}}i)"
        else:
            return f"{re_part}"
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
            
def main():
    test = "-1"
    res = str_to_latex_frac(test)

    print("original:", test)
    print("result:", res)


if __name__ == "__main__":
    main()