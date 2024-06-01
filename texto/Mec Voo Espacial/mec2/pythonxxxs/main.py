import numpy as np
from scipy.integrate import solve_ivp

from dinamica_foguete import *
from RVELPOLAR2RVELRET import *
from long_ECEF2ECI import *

global RE, WE, MUT, J2, J3, J4, G, LC, DT, R, FC, ML
global MS, M0, MP, TI, TQ, TS, ISP, H0, L_RAIL, TG, AGSO, TQ3, TQ31, TQ32, TS3, VGSO, MP3
global HAS_APOGEE, PHII_SIGN, SOL

## PROPULSÃO
ISP = [251.0, 271.0, 315.0]         # s     IMPULSO ESPECÍFICO DOS ESTÁGIOS
MP  = [5.5262E4, 11058.0, 0, 0]     # kg    MASSA DE PROPELENTE [1º, 2º] ESTÁGIO
MP3 = 243.6                         # kg    MASSA DE PROPELENTE 3º ESTÁGIO
TQ1 = 62.0                          # s     TEMPO DE QUEIMA 1º ESTÁGIO
TQ2 = 64.62                         # s     TEMPO DE QUEIMA 2º ESTÁGIO
TQ3 = 301.0                         # s     TEMPO DE QUEIMA 3º ESTÁGIO

## ESTRUTURA E CARGA UTIL
MS  = [7750.0, 1367.0, 64.7544]     # kg    MASSA ESTRUTURAL DOS ESTÁGIOS
ML  = 13.0                          # kg    MASSA DA CARGA ÚTIL

## AERODINÂMICA E AMBIENTE
FC  = 1.28                          # 1     MULTIPLICADOR DO ARRASTO

S1  = 4.6*5.0/3.0                   # m²    ÁREA DA SEÇÃO TRANSV. 1º ESTÁGIO
S2  = 1.50                          # m²    ÁREA DA SEÇÃO LONGITD. 2º ESTÁGIO
S3  = 1.50                          # m²    ÁREA DA SEÇÃO LONGITD. 3° ESTÁGIO
SL  = 1.50                          # m²    ÁREA DA SEÇÃO LONGITD. CARGA ÚTIL

LT  = 7.33 + 7.1 + 6.28             # m     COMPRIMENTO TOTAL
L2  =        7.1 + 6.28             # m     COMPRIMENTO SEM 1º ESTÁGIO
L3  =              6.28             # m     COMPRIMENTO SEM 2° ESTÁGIO
L4  = 1.0                           # m     COMprIMENTO DA CARGA ÚTIL

F2  = (L2/LT)*0.5 + 0.5             # 1     FATOR DE CORREÇÃO 1° ESTÁGIO
F3  = (L3/LT)*0.5 + 0.5             # 1     FATOR DE CORREÇÃO 2º ESTÁGIO
F4  = (L4/LT)*0.5 + 0.5             # 1     FATOR DE CORREÇÃO DA CARGA UTIL

# VETOR DE ÁREAS DE REFERÊNCIA PARA CÁLCULO DO ARRASTO (m²):
SR  = np.array([S1,     \
                S2*F2,  \
                S3*F3,  \
                SL*F4   ])

LC  = 1.50                          # m     COMPRIMENTO CARACTERÍSTICO
DT  = 10.0                          # K     DESVIO TEMP. ATMOSFERA PADRÃO

## PLANETA
RE  = 6378137.00                    # m     RAIO EQUATORIAL
WE  = 7.2921150E-5                  # rad/s VELOCIDADE INERCIAL DE ROTAÇÃO
G   = 9.80665                       # m/s²  ACELERAÇÃO GRAVITACIONAL NÍV. MAR
MUT = 3.986004418E14                # m³/s² PARÂMETRO GRAVITACIONAL

J2  = 0.00108263                    # 1     PARÂMETROS DE JEFFERY
J3  = -0.00000254                   # 1     "
J4  = -0.00000161                   # 1     "
TG  = 0.00                          # s     TEMPO DE LONGITUDE CELESTE NULA

## CONDIÇÕES INICIAIS
H0  = 0.00                          # m     ALTITUDE DA BASE DE LANÇAMENTO
DELTA0 = np.radians(-2.3267844)     # rad   LATITUDE DA BASE
LON0 = np.radians(-44.4111042)      # rad   LONGITUDE INICIAL
L_RAIL = LT                         # m     COMPRIMENTO DO TRILHO

## PARÂMETROS ORBITAIS
INGSO = np.radians(5.0)             # rad   INCLINAÇÃO
AGSO = 42.164140e6                  # m     RAIO ORBITAL
VGSO = np.sqrt(MUT / AGSO)          # m/s   VELOCIDADE ORBITAL

## PARÂMETROS PROPULSIVOS E TEMPORAIS CALCULADOS
TS1 = 2.0                           # s     TEMPO DE ESPERA ENTRE 1° e 2° ESTG.
TS2 = 2.0                           # s     TEMPO ENTRE 2º E 3º ESTÁGIO
TS3 = 2.0                           # s     TEMPO ENTRE 3° ESTÁGIO E CARGA ÚTIL

TEQ2 = 5.0                          # s     TEMPO DE ESPERA PARA IGNIÇÃO 2° EST.
TEQ3 = 900.0                        # s     TEMPO PARA IGNIÇÃO 3º ESTÁGIO

TQ31 = 251.5                        # s     DURAÇÃO DE QUEIMA 1ª IGN. 3º ESTÁGIO
TQ32 = 49.5                         # s     DURAÇÃO DE QUEIMA 2ª IGN. 3º ESTÁGIO

MP31 = MP3 * TQ31 / TQ3             # kg    MASSA PROPLNT. 3º ESTÁGIO 1ª QUEIMA
MP32 = MP3 * TQ32 / TQ3             # kg    MASSA PROPLNT. 3º ESTÁGIO 2ª QUEIMA

TI = np.array([0, 0, 0, 0])         # s     TEMPOS DE IGNIÇÃO
TQ = np.array([0, 0, 0, 0])         # s     TEMPOS DE FIM DE QUEIMA
TS = np.array([0, 0, 0])            # s     TEMPOS DE SEPARAÇÃO

TI[0] = 0.0
TQ[0] = TI[0] + TQ1
TS[0] = TQ[0] + TS1

TI[1] = TS[0] + TEQ2
TQ[1] = TI[1] + TQ2
TS[1] = TQ[1] + TS2

TI[2] = TS[1] + TEQ3
TQ[2] = TI[2] + TQ31
TS[2] = 1.0E6 + TQ32 + TS3

TI[3] = 1.0E6
TQ[3] = TI[3] + TQ32

MP[2] = MP31
MP[3] = MP32

HAS_APOGEE = False
PHII_SIGN = 0

## PARÂMETROS INICIAIS CALCULADOS
M0 = np.sum(MP) + np.sum(MS) + ML   # kg    MASSA INICIAL DO FOGUETE
R0 = RE + H0                        # m     DISTÂNCIA RADIAL INICIAL

## PREVISÃO DE DELTA-V PELA EQUAÇÃO DE FOGUETE (DV, m/s)
MPX = np.array([MP[0], MP[1], MP3]) # kg    VETOR PARA CÁLCULO DE SIGMA
SIGMA = np.divide(MS, (MS + MPX))   # 1     RAZÃO DE MASSAS

M01 = M0                            # kg    MASSA TOTAL NA DECOLAGEM
M02 = MS[1] + MPX[1] + MS[2]    \
      + MPX[2] + ML                 # kg    MASSA TOTAL NA IGN. DO 2º ESTÁGIO

M03 = MS[2] + MPX[2] + ML           # kg    MASSA TOTAL NA IGN. DO 3º ESTÁGIO


LAMB = np.array([   M02/M01,    \
                    M03/M02,    \
                    ML/M03      ])  # 1     RAZÕES DE CARGA ÚTIL

LAMBL = np.prod(LAMB)               # 1     RAZÃO DE CARGA ÚTIL TOTAL

VE = np.multiply(G, ISP)            # m/s   VELOCIDADE DE EXAUSTÃO

DV = -np.sum(                                               \
            np.multiply(VE,                                 \
                        np.log(SIGMA +                      \
                               np.multiply(1.0-SIGMA,LAMB)  \
                               )                            \
                        )                                   \
            )

################################## SIMULAÇÃO ###################################
def main():
    SIMULA = True

    while(SIMULA == True):
        
        ## ENTRADA DE PARÂMETROS PELO USUÁRIO
        TF      = input("TEMPO DE SIMULAÇÃO (s): ")
        V0      = input("VELOCIDADE RELATIVA INICIAL (m/s): ")
        PHI0    = input("ÂNGULO DE ELEVAÇÃO INICIAL (deg): ")
        PHI0    = np.radians(float(PHI0))
        
        ## AVALIAÇÃO DE VIABILIDADE
        Y = np.cos(INGSO) / np.cos(DELTA0)
        if(abs(Y) > 1):
            print("INCLINAÇÃO INATINGÍVEL PELAS COND. INICIAIS, USANDO A MENOR POSSÍVEL")
            Y = sign(Y)
        
        ## CÁLCULO DE PARÂMETROS
        AI_F = np.arcsin(Y)     # AZIMUTE FINAL DE VELOCIDADE INERCIAL
        RPGTO = RE + 250.0E3    # APOGEU DA ÓRBITA DE TRANSFERÊNCIA COM 250km ALTTD.
        AGTO = (AGSO + RPGTO)/2.0 # SEMIEIXO MAIOR ÓRBITA DE TRANSFERÊNCIA
        VIGTO = np.sqrt(MUT * (2.0/RPGTO - 1.0/AGTO)) # VEL. INERCL. ÓRBITA DE TRSF.
        
        # AZIMUTE RELATIVO:
        A0 = np.arctan(np.tan(AI_F)-(RPGTO*WE*np.cos(DELTA0))/(VIGTO*np.cos(AI_F)))

        # VETOR ESTADO INICIAL:
        X0 = np.array([V0, A0, PHI0, R0, DELTA0, LON0])
        
        # OPÇÕES DO SOLVER
        SLVR_OPT = {'rtol': 1e-8, 'atol': 1e-10, 'max_step': 0.5}
        SLVR_TSPAN = (0, TF)
        SOL = solve_ivp(dinamica_foguete, SLVR_TSPAN, X0, **SLVR_OPT)

        # EXPORTA RESULTADOS
        np.savetxt("OUT_V.DAT",         SOL.y[0], '%s')
        np.savetxt("OUT_A.DAT",         SOL.y[1], '%s')
        np.savetxt("OUT_PHI.DAT",       SOL.y[2], '%s')
        np.savetxt("OUT_R.DAT",         SOL.y[3], '%s')
        np.savetxt("OUT_DELTA0.DAT",    SOL.y[4], '%s')
        np.savetxt("OUT_LONG0.DAT",     SOL.y[5], '%s')
        np.savetxt("OUT_TIME.DAT",      SOL.t, '%s')

        ## PÓS-PROCESSAMENTO
        N = len(SOL.t)      # NÚMERO DE TIMESTEPS
        V = np.zeros(N)     # MAGNTD. VELOCIDADE
        A = np.zeros(N)     # AZIMUTE VELOCIDADE
        PHI = np.zeros(N)   # ELEVAÇÃO VELOCIDADE
        
        R = np.zeros(N)     # RAIO DE ÓRBITA
        H = np.zeros(N)     # ALTITUDE REF. PLANETA
        DELTA = np.zeros(N) # LONGITUDE REF. PLANETA
        LON = np.zeros(N)   # LONGITUDE REF. PLANETA
        
        M = np.zeros(N)     # MASSA DO FOGUETE
        
        FT = np.zeros(N)    # FORÇA PROPULSIVA
        MU = np.zeros(N)    # ÂNGULO PROPULSIVO MU
        EPSL = np.zeros(N)  # ÂNGULO PROPULSIVO EPSILON
        
        D = np.zeros(N)     # FORÇA DE ARRASTO
        Q = np.zeros(N)     # PRESSÃO DINÂMICA
        MA = np.zeros(N)    # NÚMERO DE MACH
        KN = np.zeros(N)    # NÚMERO DE KNUDSEN
        T = np.zeros(N)     # TEMPERATURA
        RHO = np.zeros(N)   # MASSA ESPECÍFICA ATMOSFÉRICA
        VI = np.zeros(N)    # MAGNTD. VELOCIDADE INERCIAL
        AI = np.zeros(N)    # AZIMUTE VELOCIDADE INERCIAL
        PHII = np.zeros(N)  # ELEVAÇÃO VELOCIDADE INERCIAL
        
        LONGC = np.zeros(N) # LONGITUDE CELESTE
        EE = np.zeros(N)    # ENERGIA ESPECÍFICA
        A = np.zeros(N)     # SEMIEIXO MAIOR ORBITAL
        E = np.zeros(N)     # EXCENTRICIDADE ORBITAL
        TAU = np.zeros(N)   # TEMPO DE PERIGEU
        OM = np.zeros(N)    # ASCENÇÃO RETA DO NODO ASCENDENTE
        IN = np.zeros(N)    # INCLINAÇÃO ORBITAL
        OM = np.zeros(N)    # ARGUMENTO DE PERIGEU
        R0b = np.zeros(N)   # POSIÇÃO NO REFERENCIAL ECI
        
        for i in range(N):
            V[i]        = float(SOL.y[0][i])
            A[i]        = float(SOL.y[1][i])
            PHI[i]      = float(SOL.y[2][i])
            R[i]        = float(SOL.y[3][i])
            DELTA[i]    = float(SOL.y[4][i])
            LON[i]      = float(SOL.y[5][i])
            
            H[i]        = R[i] - RE
            
            TI = float(SOL.t[i])
            XI= [V[i], A[i], PHI[i], R[i], DELTA[i], LON[i]]
            
            FT[i], M[i], MU[i], EPSL[i] = propulsao_N_estagios(TI, XI)
            T[i], _, _, RHO[i], _, MA[i], _, _, KN[i], _, _, R[i] = atm_padrao(H[i], V[i], LC, DT)
            D[i], _, _ = aerodinamica_N_estagios(TI, V[i], H[i], MA[i], KN[i], T[i], RHO[i], R[i])
            Q[i] = 0.5 * RHO[i] * V[i]**2
            VI[i], PHII[i], AI[i] = vRel2vIne(V[i], PHII[i], A[i], WE, R[i], DELTA[i])
            LONGC[i] = long_ECEF2ECI(TI, LON[i], WE, TG)
            EE[i] = (VI[i]**2) / 2 - MUT/R[i]
            RC0, VC0 = RvelPolar2RvelRet(VI[i], AI[i], PHII[i], R[i], DELTA[i], LONGC[i])
            #print(RC0)
            #R0[i] = RC0.T

        
        SIMULA = input('Deseja simular novamente (s/n): ').lower().strip() == 's'

if __name__ == '__main__':
    main()
