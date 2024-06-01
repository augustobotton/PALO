import numpy as np

from main import *

def propulsao_N_estagios(t, X):
    global TI, TQ, TS, ISP, MP, MS, M0, G, RE

    N = len(TS)
    
    if N == 1:
        ft, m = propulsor_1_estagio(t, TI, TQ, TS, ISP, MP, MS, M0, G)
    elif N == 2:
        ft, m = propulsor_2_estagios(t, TI, TQ, TS, ISP, MP, MS, M0, G)
    elif N == 3:
        ft, m = propulsor_3_estagios(t, TI, TQ, TS, ISP, MP, MS, M0, G)

    V = X[0]
    A = X[1]
    PHI = X[2]
    R = X[3]
    DELTA = X[4]
    
    H = R - RE
    EPSL = 0.0
    MU = 0.0

    if(H > 200.0E3):
        _, PHII, AI = Vrel2Vine(V, PHI, A, WE, R, DELTA)
        MU=anp.sin(np.cos(A)*np.cos(PHII)*np.sin(AI)-np.sin(A)*np.cos(PHII)*np.cos(AI))
        
        EPSL=atan2(-np.cos(PHI)*np.sin(PHII)+np.sin(PHI)*np.sin(A)*np.cos(PHII)*np.sin(AI)+			\
                 np.sin(PHI)*np.cos(A)*np.cos(PHII)*np.cos(AI),np.sin(PHI)*np.sin(PHII)+			\
                 np.cos(PHI)*np.sin(A)*np.cos(PHII)*np.sin(AI)+np.cos(PHI)*np.cos(A)*np.cos(PHII)*np.cos(AI))
        
    return ft, m, MU, EPSL

def propulsor_1_estagio(t, TI, TQ, TS, ISP, MP, MS, M0, G):
    if t <= TI[0]:
        m = M0
        ft = 0
    elif t <= TQ[0]:
        md = -MP[0] / (TQ[0] - TI[0])
        m = M0 + md * (t - TI[0])
        ft = -G * ISP[0] * md
    elif t <= TS[0]:
        m = M0 - MP[0]
        ft = 0
    else:
        m = M0 - MP[0] - MS[0]
        ft = 0
    
    return ft, m

def propulsor_2_estagios(t, TI, TQ, TS, ISP, MP, MS, M0, G):
    if t <= TI[0]:
        m = M0
        ft = 0
    elif t <= TQ[0]:
        md = -MP[0] / (TQ[0] - TI[0])
        m = M0 + md * (t - TI[0])
        ft = -G * ISP[0] * md
    elif t <= TS[0]:
        m = M0 - MP[0]
        ft = 0
    elif t <= TI[1]:
        m = M0 - MP[0] - MS[0]
        ft = 0
    elif t <= TQ[1]:
        md = -MP[1] / (TQ[1] - TI[1])
        M02 = M0 - MP[0] - MS[0]
        m = M02 + md * (t - TI[1])
        ft = -G * ISP[1] * md
    elif t <= TS[1]:
        m = M0 - MP[0] - MS[0] - MP[1]
        ft = 0
    else:
        m = M0 - MP[0] - MS[0] - MP[1] - MS[1]
        ft = 0
    
    return ft, m


def propulsor_3_estagios(t, TI, TQ, TS, ISP, MP, MS, M0, G):
    if t <= TI[0]:
        m = M0
        ft = 0
    elif t <= TQ[0]:
        md = -MP[0] / (TQ[0] - TI[0])
        m = M0 + md * (t - TI[0])
        ft = -G * ISP[0] * md
    elif t <= TS[0]:
        m = M0 - MP[0]
        ft = 0
    elif t <= TI[1]:
        m = M0 - MP[0] - MS[0]
        ft = 0
    elif t <= TQ[1]:
        md = -MP[1] / (TQ[1] - TI[1])
        M02 = M0 - MP[0] - MS[0]
        m = M02 + md * (t - TI[1])
        ft = -G * ISP[1] * md
    elif t <= TS[1]:
        m = M0 - MP[0] - MS[0] - MP[1]
        ft = 0
    elif t <= TI[2]:
        m = M0 - MP[0] - MS[0] - MP[1] - MS[1]
        ft = 0
    elif t <= TQ[2]:
        md = -MP[2] / (TQ[2] - TI[2])
        M03 = M0 - MP[0] - MS[0] - MP[1] - MS[1]
        m = M03 + md * (t - TI[2])
        ft = -G * ISP[2] * md
    elif t <= TS[2]:
        m = M0 - MP[0] - MS[0] - MP[1] - MS[1] - MP[2]
        ft = 0
    else:
        m = M0 - MP[0] - MS[0] - MP[1] - MS[1] - MP[2] - MS[2]
        ft = 0

    return ft, m
