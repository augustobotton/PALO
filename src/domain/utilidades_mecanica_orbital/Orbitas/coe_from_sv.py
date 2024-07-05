import numpy as np


def coe_from_sv(R, V, mu):
    """
    Calcula os elementos orbitais clássicos (coe) a partir do vetor de estado (R, V).
    """
    eps = 1.e-10

    r = np.linalg.norm(R)
    v = np.linalg.norm(V)

    vr = np.dot(R, V) / r

    H = np.cross(R, V)
    h = np.linalg.norm(H)

    incl = np.arccos(H[2] / h)

    K = np.array([0, 0, 1])
    N = np.cross(K, H)
    n = np.linalg.norm(N)

    if n != 0:
        RA = np.arccos(N[0] / n)
        if N[1] < 0:
            RA = 2 * np.pi - RA
    else:
        RA = 0

    E = (1 / mu) * ((v ** 2 - (mu / r)) * R - (r * vr * V))
    e = np.linalg.norm(E)

    if n != 0:
        if e > eps:
            w = np.arccos(np.dot(N, E) / (n * e))
            if E[2] < 0:
                w = 2 * np.pi - w
        else:
            w = 0
    else:
        w = 0

    if e > eps:
        TA = np.arccos(np.dot(E, R) / (e * r))
        if vr < 0:
            TA = 2 * np.pi - TA
    else:
        cp = np.cross(N, R)
        if cp[2] >= 0:
            TA = np.arccos(np.dot(N, R) / (n * r))
        else:
            TA = 2 * np.pi - np.arccos(np.dot(N, R) / (n * r))

    a = h ** 2 / mu / (1 - e ** 2)

    return a, e, incl, RA, w, TA, h
