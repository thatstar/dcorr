"""
Lebedev-Laikov quadrature grids
"""

import numpy as np


def genOh_a00(v):
    "(0, 0, a) etc. (6 points)"
    a = 1.0
    r = [( a,  0,  0, v),
         (-a,  0,  0, v),
         ( 0,  a,  0, v),
         ( 0, -a,  0, v),
         ( 0,  0,  a, v),
         ( 0,  0, -a, v)]

    return r


def genOh_aa0(v):
    "(0, a, a) etc, a = 1/np.sqrt(2) (12 points)"
    a = 1/np.sqrt(2)
    r = [( 0,  a,  a, v),
         ( 0, -a,  a, v),
         ( 0,  a, -a, v),
         ( 0, -a, -a, v),
         ( a,  0,  a, v),
         (-a,  0,  a, v),
         ( a,  0, -a, v),
         (-a,  0, -a, v),
         ( a,  a,  0, v),
         (-a,  a,  0, v),
         ( a, -a,  0, v),
         (-a, -a,  0, v)]

    return r


def genOh_aaa(v):
    "(a, a, a) etc, a = 1/np.sqrt(3) (8 points)"
    a = np.sqrt(1./3.)
    r = [( a,  a,  a, v),
         (-a,  a,  a, v),
         ( a, -a,  a, v),
         (-a, -a,  a, v),
         ( a,  a, -a, v),
         (-a,  a, -a, v),
         ( a, -a, -a, v),
         (-a, -a, -a, v)]

    return r


def genOh_aab(v, a):
    "(a, a, b) etc, b = np.sqrt(1 - 2*a**2), a input (24 points)"
    b = np.sqrt(1.0 - 2.0*a*a)
    r = [( a,  a,  b, v),
         (-a,  a,  b, v),
         ( a, -a,  b, v),
         (-a, -a,  b, v),
         ( a,  a, -b, v),
         (-a,  a, -b, v),
         ( a, -a, -b, v),
         (-a, -a, -b, v),
         ( a,  b,  a, v),
         (-a,  b,  a, v),
         ( a, -b,  a, v),
         (-a, -b,  a, v),
         ( a,  b, -a, v),
         (-a,  b, -a, v),
         ( a, -b, -a, v),
         (-a, -b, -a, v),
         ( b,  a,  a, v),
         (-b,  a,  a, v),
         ( b, -a,  a, v),
         (-b, -a,  a, v),
         ( b,  a, -a, v),
         (-b,  a, -a, v),
         ( b, -a, -a, v),
         (-b, -a, -a, v)]

    return r


def genOh_ab0(v,a):
    "(a, b, 0) etc, b = np.sqrt(1 - a**2), a input (24 points)"
    b = np.sqrt(1.0 - a*a)
    r = [( a,  b,  0, v),
         (-a,  b,  0, v),
         ( a, -b,  0, v),
         (-a, -b,  0, v),
         ( b,  a,  0, v),
         (-b,  a,  0, v),
         ( b, -a,  0, v),
         (-b, -a,  0, v),
         ( a,  0,  b, v),
         (-a,  0,  b, v),
         ( a,  0, -b, v),
         (-a,  0, -b, v),
         ( b,  0,  a, v),
         (-b,  0,  a, v),
         ( b,  0, -a, v),
         (-b,  0, -a, v),
         ( 0,  a,  b, v),
         ( 0, -a,  b, v),
         ( 0,  a, -b, v),
         ( 0, -a, -b, v),
         ( 0,  b,  a, v),
         ( 0, -b,  a, v),
         ( 0,  b, -a, v),
         ( 0, -b, -a, v)]

    return r


def genOh_abc(v,a,b):
    "(a, b, c) etc, c = np.sqrt(1 - a**2 - b**2), a, b input (48 points)"
    c = np.sqrt(1.0 - a*a - b*b)
    r = [( a,  b,  c, v),
         (-a,  b,  c, v),
         ( a, -b,  c, v),
         (-a, -b,  c, v),
         ( a,  b, -c, v),
         (-a,  b, -c, v),
         ( a, -b, -c, v),
         (-a, -b, -c, v),
         ( a,  c,  b, v),
         (-a,  c,  b, v),
         ( a, -c,  b, v),
         (-a, -c,  b, v),
         ( a,  c, -b, v),
         (-a,  c, -b, v),
         ( a, -c, -b, v),
         (-a, -c, -b, v),
         ( b,  a,  c, v),
         (-b,  a,  c, v),
         ( b, -a,  c, v),
         (-b, -a,  c, v),
         ( b,  a, -c, v),
         (-b,  a, -c, v),
         ( b, -a, -c, v),
         (-b, -a, -c, v),
         ( b,  c,  a, v),
         (-b,  c,  a, v),
         ( b, -c,  a, v),
         (-b, -c,  a, v),
         ( b,  c, -a, v),
         (-b,  c, -a, v),
         ( b, -c, -a, v),
         (-b, -c, -a, v),
         ( c,  a,  b, v),
         (-c,  a,  b, v),
         ( c, -a,  b, v),
         (-c, -a,  b, v),
         ( c,  a, -b, v),
         (-c,  a, -b, v),
         ( c, -a, -b, v),
         (-c, -a, -b, v),
         ( c,  b,  a, v),
         (-c,  b,  a, v),
         ( c, -b,  a, v),
         (-c, -b,  a, v),
         ( c,  b, -a, v),
         (-c,  b, -a, v),
         ( c, -b, -a, v),
         (-c, -b, -a, v)]

    return r


def leb6():
    leb = genOh_a00(0.1666666666666667)

    return leb


def leb14():
    leb = genOh_a00(0.06666666666666667) \
        + genOh_aaa(0.07500000000000000)

    return leb


def leb26():
    leb = genOh_a00(0.04761904761904762) \
        + genOh_aa0(0.03809523809523810) \
        + genOh_aaa(0.03214285714285714)

    return leb


def leb38():
    leb = genOh_a00(0.009523809523809524) \
        + genOh_aaa(0.032142857142857140) \
        + genOh_ab0(0.028571428571428570, 0.4597008433809831)

    return leb


def leb50():
    leb = genOh_a00(0.01269841269841270) \
        + genOh_aa0(0.02257495590828924) \
        + genOh_aaa(0.02109375000000000) \
        + genOh_aab(0.02017333553791887, 0.3015113445777636)

    return leb


def leb74():
    leb = genOh_a00( 0.5130671797338464E-3) \
        + genOh_aa0( 0.1660406956574204E-1) \
        + genOh_aaa(-0.2958603896103896E-1) \
        + genOh_aab( 0.2657620708215946E-1, 0.4803844614152614E+0) \
        + genOh_ab0( 0.1652217099371571E-1, 0.3207726489807764E+0)

    return leb


def leb86():
    leb = genOh_a00(0.1154401154401154E-1) \
        + genOh_aaa(0.1194390908585628E-1) \
        + genOh_aab(0.1111055571060340E-1, 0.3696028464541502E+0) \
        + genOh_aab(0.1187650129453714E-1, 0.6943540066026664E+0) \
        + genOh_ab0(0.1181230374690448E-1, 0.3742430390903412E+0)

    return leb


def leb110():
    leb = genOh_a00(0.3828270494937162E-2) \
        + genOh_aaa(0.9793737512487512E-2) \
        + genOh_aab(0.8211737283191111E-2, 0.1851156353447362E+0) \
        + genOh_aab(0.9942814891178103E-2, 0.6904210483822922E+0) \
        + genOh_aab(0.9595471336070963E-2, 0.3956894730559419E+0) \
        + genOh_ab0(0.9694996361663028E-2, 0.4783690288121502E+0)

    return leb


def leb146():
    leb = genOh_a00(0.5996313688621381E-3) \
        + genOh_aa0(0.7372999718620756E-2) \
        + genOh_aaa(0.7210515360144488E-2) \
        + genOh_aab(0.7116355493117555E-2,0.6764410400114264E+0) \
        + genOh_aab(0.6753829486314477E-2,0.4174961227965453E+0) \
        + genOh_aab(0.7574394159054034E-2,0.1574676672039082E+0) \
        + genOh_abc(0.6991087353303262E-2,0.1403553811713183E+0, 0.4493328323269557E+0)

    return leb


def leb170():
    leb = genOh_a00(0.5544842902037365E-2) \
        + genOh_aa0(0.6071332770670752E-2) \
        + genOh_aaa(0.6383674773515093E-2) \
        + genOh_aab(0.5183387587747790E-2, 0.2551252621114134E+0) \
        + genOh_aab(0.6317929009813725E-2, 0.6743601460362766E+0) \
        + genOh_aab(0.6201670006589077E-2, 0.4318910696719410E+0) \
        + genOh_ab0(0.5477143385137348E-2, 0.2613931360335988E+0) \
        + genOh_abc(0.5968383987681156E-2, 0.4990453161796037E+0, 0.1446630744325115E+0)

    return leb


def leb194():
    leb = genOh_a00(0.1782340447244611E-2) \
        + genOh_aa0(0.5716905949977102E-2) \
        + genOh_aaa(0.5573383178848738E-2) \
        + genOh_aab(0.5608704082587997E-2, 0.6712973442695226E+0) \
        + genOh_aab(0.5158237711805383E-2, 0.2892465627575439E+0) \
        + genOh_aab(0.5518771467273614E-2, 0.4446933178717437E+0) \
        + genOh_aab(0.4106777028169394E-2, 0.1299335447650067E+0) \
        + genOh_ab0(0.5051846064614808E-2, 0.3457702197611283E+0) \
        + genOh_abc(0.5530248916233094E-2, 0.1590417105383530E+0, 0.8360360154824589E+0)

    return leb


def lebedev_grid(n):
    _funcs_list = {
          6: leb6,
         14: leb14,
         26: leb26,
         38: leb38,
         50: leb50,
         74: leb74,
         86: leb86,
        110: leb110,
        146: leb146,
        170: leb170,
        194: leb194
    }

    return np.array(_funcs_list[n]())


if __name__ == "__main__":
    import sys
    r = float(sys.argv[1])
    n = int(sys.argv[2])
    s = lebedev_grid(n)
    for i in range(s.shape[0]):
        print(s[i][0]*r, s[i][1]*r, s[i][2]*r, s[i][3])
