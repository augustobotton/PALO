from VREL2VINE import *
import math
import numpy as np

from main import *

def parametros_manobra_adquire_gso(t, m, X):
    #   ENTRADA
    #   VARIÁVEL    UNIDADE DESCRIÇÃO
    #   t           s       Tempo
    #   m           kg      Massa do foguete
    #   X[0]        m/s     Módulo vetor velocidade
    #   X[1]        rad     Ângulo de azimute vetor velocidade
    #   X[2]        rad     Ângulo de elevação vetor velocidade
    #   X[3]        m       Distância radial até o centro de massa do planeta
    #   X[4]        rad     Latitude ref. plano equatorial do planeta
    #   X[5]        rad     Longitude planetária
    
    global WE, AGSO, PHII_SIGN, HAS_APOGEE, TI, TQ, TS, TQ3, TQ31, TQ32, TS3, VGSO, ISP, G, MP3, WE
    
    V       = X[0]
    A       = X[1]
    PHI     = X[2]
    R       = X[3]
    DELTA   = X[4]
    LON     = X[5]

    VI, PHII, _ = vRel2vIne(V, PHI, A, WE, R, DELTA)

    if( (R > 0.9*AGSO)          and     \
        (HAS_APOGEE == False)   and     \
        (sign(PHII) != PHII_SIGN)       ):
        
        HAS_APOGEE = True   # Encontrou o apogeu
        TI[3] = t           # Iniciar segunda ignição
        
        # Cálculo do tempo de queima para o dV desejado
        DVGSO = VGSO - VI
        MP32 = (M*math.exp(DVGSO / (ISP[2] * G)) - m)   \
               / math.exp(DVGSO / (ISP[2] * G))
        
        TQ32 = TQ3 * MP32 / MP3 # Tempo de queima
        
        if(TQ31 + TQ32 > TQ3):
            TQ32 = TQ3 - TQ31
        
        TQ[3] = TI[3] + TQ32
        TS[3] = TQ[3] + TS3
        
    PHII_SIGN = np.sign(PHII)
