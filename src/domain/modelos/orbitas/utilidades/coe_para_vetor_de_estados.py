import numpy as np


def coe_para_vetor_de_estados(coe, mu):
    """
    This function computes the state vector (r, v) from the classical orbital elements (coe).
    """

    e = coe[1]
    RA = coe[3]
    incl = coe[2]
    w = coe[4]
    TA = coe[5]
    h = coe[6]

    rp = (h ** 2 / mu) * (1 / (1 + e * np.cos(TA))) * np.array([np.cos(TA), np.sin(TA), 0])
    vp = (mu / h) * np.array([-np.sin(TA), e + np.cos(TA), 0])

    R3_RA = np.array([
        [np.cos(RA), np.sin(RA), 0],
        [-np.sin(RA), np.cos(RA), 0],
        [0, 0, 1]
    ])

    R1_i = np.array([
        [1, 0, 0],
        [0, np.cos(incl), np.sin(incl)],
        [0, -np.sin(incl), np.cos(incl)]
    ])

    R3_w = np.array([
        [np.cos(w), np.sin(w), 0],
        [-np.sin(w), np.cos(w), 0],
        [0, 0, 1]
    ])

    Q_pX = (R3_w @ R1_i @ R3_RA).T

    r = Q_pX @ rp
    v = Q_pX @ vp

    return r, v
