import numpy as np
from modelo_aerodinamico import *
from main import *

def aerodinamica_N_estagios(t, V, h, M, Kn, T, rho, R):
    global TS, SR, FC

    CD = modelo_aerodinamico(V, h, M, Kn, T, R)
    CD = FC * CD

    N = len(TS)
    if N == 1:
        S = area1estagio(t, TS, SR)
    elif N == 2:
        S = area2estagios(t, TS, SR)
    elif N == 3:
        S = area3estagios(t, TS, SR)

    D = 0.5 * rho * V**2 * S * CD
    fy = 0
    L = 0

    return D, fy, L

def area1estagio(t, TS, SR):
    if t <= TS[0]:
        S = SR[0]
    else:
        S = SR[1]

    return S

def area2estagios(t, TS, SR):
    if t <= TS[0]:
        S = SR[0]
    elif t <= TS[1]:
        S = SR[1]
    else:
        S = SR[2]

    return S

def area3estagios(t, TS, SR):
    if t <= TS[0]:
        S = SR[0]
    elif t <= TS[1]:
        S = SR[1]
    elif t <= TS[2]:
        S = SR[2]
    else:
        S = SR[3]

    return S