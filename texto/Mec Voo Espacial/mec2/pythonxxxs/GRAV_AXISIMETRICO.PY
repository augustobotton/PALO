import numpy as np
from main import *

def grav_axisimetrico(R, DELTA):
    global RE, MUT, J2, J3, J4
    #   ENTRADAS
    #   VARIÁVEL    UNIDADE DESCRIÇÃO
    #   R           m       Distância radial ao centro de massa do planeta
    #   DELTA       rad     Latitude centrada no planeta
    #
    #   SAÍDAS
    #   VARIÁVEL    UNIDADE DESCRIÇÃO
    #   GC          m/s²    Componente centrípeta gravitacional
    #   GDEL        m/s²    Componente latitudinal gravitacional

    PHI = np.pi / 2.0 - DELTA # Co-latitude
    
    GR = -(MUT) / R ** 2 - (1 / R) * MUT * (-((J2 * RE ** 2 * (-1 + 3 * np.cos(PHI) ** 2)) / R ** 3) - (1.5 * J3 * RE ** 3 * (-3 * np.cos(PHI) + 5 * np.cos(PHI) ** 3)) / R ** 4 - (J4 * RE ** 4 * (3 - 30 * np.cos(PHI) ** 2 + 35 * np.cos(PHI) ** 4)) / (2 * R ** 5)) + (MUT / R ** 2) * ((0.5 * J2 * RE ** 2 * (-1 + 3 * np.cos(PHI) ** 2)) / R ** 2 + (0.5 * J3 * RE ** 3 * (-3 * np.cos(PHI) + 5 * np.cos(PHI) ** 3)) / R ** 3 + (J4 * RE ** 4 * (3 - 30 * np.cos(PHI) ** 2 + 35 * np.cos(PHI) ** 4)) / (8 * R ** 4))
    
    GPHI = -(1 / (R ** 2)) * MUT * (-((3 * J2 * RE ** 2 * np.cos(PHI) * np.sin(PHI)) / R ** 2) + (0.5 * J3 * RE ** 3 * (3 * np.sin(PHI) - 15 * np.cos(PHI) ** 2 * np.sin(PHI))) / R ** 3 + (J4 * RE ** 4 * (60 * np.cos(PHI) * np.sin(PHI) - 140 * np.cos(PHI) ** 3 * np.sin(PHI))) / (8 * R ** 4))

    GC = -GR
    GDEL = -GPHI
    
    return GC, GDEL
