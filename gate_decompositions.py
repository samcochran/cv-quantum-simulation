import numpy as np
from math import factorial
from scipy.special import binom
import strawberryfields as sf
from strawberryfields.ops import Pgate, Vgate, CZgate, CXgate, Fouriergate
# import warnings
# warnings.filterwarnings(action='ignore', module='strawberryfields')
np.random.seed(42)
sf.hbar = 1


def quadratic2(C, t, j, k, q):
    """
    Quantum circuit for the 2-mode quadratic term e^(-1j C t P_j X_k) (just a wrapper for the CXgate)
    Parameters:
        C: constant coefficient
        t: evolution time
        j, k: qumode indices
        q: program qumodes
    """
    CXgate(t*C) | (q[k], q[j])


def cubic2(C, t, j, k, q):
    """
    Quantum circuit for the 2-mode cubic term e^(1j C t P_j X_k^2)
    Parameters:
        C: positive coefficient
        t: evolution time
        j, k: qumode indices
        q: program qumodes
    """
    alpha = np.sqrt(C/3)
    # e^(i 3alpha^3t/4 X_k^3)
    Vgate(9*alpha**3*t/4) | q[k]
    # e^(-i t P_j^3)
    Fouriergate() | q[j]
    Vgate(3*t) | q[j]
    Fouriergate().H | q[j]
    # e^(i alpha X_k X_j)
    CZgate(alpha) | (q[k], q[j])
    # e^(i t P_j^3)
    Fouriergate() | q[j]
    Vgate(-3*t) | q[j]
    Fouriergate().H | q[j]
    # e^(-i alpha X_k X_j)
    CZgate(-alpha) | (q[k], q[j])
    # e^(-i t P_j^3)
    Fouriergate() | q[j]
    Vgate(3*t) | q[j]
    Fouriergate().H | q[j]
    # e^(-i alpha X_k X_j)
    CZgate(-alpha) | (q[k], q[j])
    # e^(i t P_j^3)
    Fouriergate() | q[j]
    Vgate(-3*t) | q[j]
    Fouriergate().H | q[j]
    # e^(i alpha X_k X_j)
    CZgate(alpha) | (q[k], q[j])


def cubic2X(C, t, j, k, q):
    """
    Quantum circuit for the 2-mode cubic term e^(-1j C t X_j X_k^2)
    Parameters:
        C: positive coefficient
        t: evolution time
        j, k: qumode indices
        q: program qumodes
    """
    Fouriergate() | q[j]
    cubic2(C, -t, j, k, q)
    Fouriergate().H | q[j]


def cubic3X(C, t, j, k, l, q):
    """
    Quantum circuit for the 3-mode cubic term e^(1j C t X_j X_k X_l)
    Parameters:
        C: constant coefficient
        t: evolution time
        j, k, l: qumode indices
        q: program qumodes
    """
    alpha = C*t

    def identity1(alpha, j, k):
        # e^(i alpha (X_j + X_k)^3)
        CXgate(1) | (q[k], q[j])
        Vgate(3*alpha) | q[j]
        CXgate(-1) | (q[k], q[j])

    def identity2(alpha, j, k, l):
        # e^(i alpha (X_j + X_k + X_l)^3)
        CXgate(1) | (q[l], q[j])
        identity1(alpha, j, k)
        CXgate(-1) | (q[l], q[j])

    Vgate(alpha/2) | q[l]
    Vgate(alpha/2) | q[k]
    Vgate(alpha/2) | q[j]
    identity1(-alpha/6, k, l)
    identity1(-alpha/6, j, l)
    identity1(-alpha/6, j, k)
    identity2(alpha/6, j, k, l)


def cubic3(C, t, j, k, l, q):
    """
    Quantum circuit for the 3-mode cubic term e^(-1j C t P_j X_k X_l)
    Parameters:
        C: constant coefficient
        t: evolution time
        j, k, l: qumode indices
        q: program qumodes
    """
    Fouriergate() | q[j]
    cubic3X(C, t, j, k, l, q)
    Fouriergate().H | q[j]


def quartic1(alpha, j, ancilla, q):
    """
    Quantum circuit for the single-mode quartic term e^(1j alpha X_j^4)
    Parameters:
        C: positive coefficient
        t: evolution time
        j, k, l: qumode indices
        q: program qumodes
    """
    def identity(alpha, j, ancilla):
        # e^(i alpha (X_j^2 + X_k)^2)
        cubic2(C=1, t=-1, j=ancilla, k=j, q=q)
        Pgate(2*alpha) | q[ancilla]
        cubic2(C=1, t=1, j=ancilla, k=j, q=q)

    cubic2X(C=2, t=alpha, j=ancilla, k=j, q=q)
    Pgate(-2*alpha) | q[ancilla]
    identity(alpha, j, ancilla)


def quarticX(a_list, C, t, j, k, l, m, ancilla, q):
    """
    Quantum circuit for the 2-4-mode quartic term e^(1j C t X_j X_k^a_1 X_l^a_2 X_m^a_3) with a_1+a_2+a_3=3
    Parameters:
        a_list: length 4 list of a coefficients [1, a_2, a_3, a_4]
        C: positive coefficient
        t: evolution time
        j, k, l: qumode indices
        ancilla: index of ancilla qumode
        q: program qumodes
    """
    assert a_list[0] == 1
    a = sum(a_list)
    assert a == 4

    def unitary(h):
        if h[3] != 0: CXgate(h[3]) | (q[m], q[j])
        if h[2] != 0: CXgate(h[2]) | (q[l], q[j])
        if h[1] != 0: CXgate(h[1]) | (q[k], q[j])

    def unitaryH(h):
        if h[1] != 0: CXgate(h[1]).H | (q[k], q[j])
        if h[2] != 0: CXgate(h[2]).H | (q[l], q[j])
        if h[3] != 0: CXgate(h[3]).H | (q[m], q[j])

    def coefficient(v):
        return 1/(2**(a-1)*factorial(a))*(-1)**sum(v[1:])*np.prod([binom(a_list[i], v[i]) for i in range(1, 4)])

    def multiplicand(v, h):
        unitary(h)
        quartic1(C*t*coefficient(v), j, ancilla, q)
        unitaryH(h)
    for v_1 in range(a_list[1] + 1):
        for v_2 in range(a_list[2] + 1):
            for v_3 in range(a_list[3] + 1):
                v = [0, v_1, v_2, v_3]
                h = [a_list[i]-2*v_i for i, v_i in zip(range(len(v)), v)]
                multiplicand(v, h)


def quartic(a_list, C, t, j, k, l, m, ancilla, q):
    """
    Quantum circuit for the 2-4-mode quartic term e^(-1j C t P_j X_k^a_2 X_l^a_3 X_m^a_4) with a+b+c=3
    Parameters:
        a_list: length 4 list of a coefficients [1, a_2, a_3, a_4]
        C: positive coefficient
        t: evolution time
        j, k, l: qumode indices
        ancilla: index of ancilla qumode
        q: program qumodes
    """
    Fouriergate() | q[j]
    quarticX(a_list, C, t, j, k, l, m, ancilla, q)
    Fouriergate() | q[j]
