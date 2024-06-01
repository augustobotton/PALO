import numpy as np

def RvelPolar2RvelRet(V, A, PHI, R, LAT, LON):
    #   ENTRADAS
    #   VARIÁVEL    UNIDADE DESCRIÇÃO
    #   V           m/s     Magnitude do vetor velocidade
    #   A           rad     Azimute da velocidade
    #   PHI         rad     Elevação da velocidade
    #   R           m       Distância radial
    #   LAT         rad     Latitude
    #   LON         rad     Longitude
    #
    #   SAÍDAS
    #   VARIÁVEL    UNIDADE DESCRIÇÃO
    #   R           m   (3) Vetor posição coords. retangulares ECI/ECEF
    #   V           m/s (3) Vetor velocidade coords. retangualres ECI/ECEF

    # MATRIZ DE CONVERSÃO ECI/ECEF -> LVLH
    CLH=np.array([                                                          \
        [np.cos(LAT)*np.cos(LON), np.cos(LAT)*np.sin(LON), np.sin(LAT)],    \
        [-np.sin(LON), np.cos(LON), 0.0],                                   \
        [-np.sin(LAT)*np.cos(LON), -np.sin(LAT)*np.sin(LON), np.cos(LAT)]   \
    ])

    # VETOR VELOCIDADE LVLH
    VLVLH=V*np.array([              \
        [np.sin(PHI)],              \
        [np.cos(PHI)*np.sin(A)],    \
        [np.cos(PHI)*np.cos(A)]     \
    ])
    
    # VELOCIDADE LVLH -> ECI/ECEF
    V=CLH.T.dot(VLVLH)

    # POSIÇÃO NO SISTEMA LVLH
    RLVLH = np.array([  [R],
                        [0],
                        [0] ])
    
    # POSIÇÃO LVLH -> ECI/ECEF
    R = CLH.T.dot(RLVLH)
    
    return R, V
