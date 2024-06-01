import numpy as np
from scipy.integrate import solve_ivp
from propulsao_N_estagios import *
from atm_padrao import *
from aerodinamica_N_estagios import *
from GRAV_AXISIMETRICO import *
from PARAMETROS_MANOBRA_ADQUIRE_GSO import *
from main import *

def dinamica_foguete(t, X):
    global WE, RE, LC, DT, H0, L_RAIL
    V, A, phi, r, delta, lon = X

    if V < 0:
        V = 0  # Evita velocidade negativa 

    # Função para cáLCulo da massa e tração em função do tempo
    ft, m, mu, epsl = propulsao_N_estagios(t, X)

    # Função para cáLCulo do modelo atmosférico
    h = r - RE  # Altitude
    T, _, _, rho, _, M, _, _, Kn, _, _, R = atm_padrao(h, V, LC, DT)

    # CáLCulo do modelo aerodinâmico
    D, fy, L = aerodinamica_N_estagios(t, V, h, M, Kn, T, rho, R)

    # CaLCulo da gravidade
    gc, gd = grav_axisimetrico(r, delta)

    # Equações de cinemática de translação
    rp = V * np.sin(phi)
    deltap = (V / r) * np.cos(phi) * np.cos(A)
    lonp = (V * np.cos(phi) * np.sin(A)) / (r * np.cos(delta))

    # Equações de dinâmica de translação
    Vp = (1 / m) * (ft * np.cos(epsl) * np.cos(mu) - D - m * gc * np.sin(phi) + m * gd * np.cos(phi) * np.cos(A) -
                    m * WE ** 2 * r * np.cos(delta) * (np.cos(phi) * np.cos(A) * np.sin(delta) - np.sin(phi) * np.cos(delta)))
    Ap = (1 / (m * V * np.cos(phi))) * (m * (V ** 2 / r) * np.cos(phi) ** 2 * np.sin(A) * np.tan(delta) + ft * np.sin(mu) + fy -
                                         m * gd * np.sin(A) +
                                         m * WE ** 2 * r * np.sin(A) * np.sin(delta) * np.cos(delta) - 2 * m * WE * V *
                                         (np.sin(phi) * np.cos(A) * np.cos(delta) - np.cos(phi) * np.sin(delta)))
    phip = (1 / (m * V)) * (m * (V ** 2 / r) * np.cos(phi) + ft * np.sin(epsl) * np.cos(mu) + L - m * gc * np.cos(phi) -
                            m * gd * np.sin(phi) * np.cos(A) +
                            m * WE ** 2 * r * np.cos(delta) * (np.sin(phi) * np.cos(A) * np.sin(delta) + np.cos(phi) * np.cos(delta)) +
                            2 * m * WE * V * np.sin(A) * np.cos(delta))

    # Saturação da altitude
    if h < 0:  # Altitude negativa não é permitida
        # Mantém as derivadas nulas
        rp = 0
        deltap = 0
        lonp = 0
        Vp = 0
        Ap = 0
        phip = 0

    # Modela o trilho de lancamento
    H = h - H0  # Altura
    if (H <= L_RAIL) and (t <= 10):  # Verifica se a altura é menor que L_RAIL nos primeiros segundos da simulação
        Ap = 0
        phip = 0  # Anula as derivadas dos ângulos de orientação da velocidade

    # Modela a navegação para o segundo disparo do motor do terceiro estágio
    parametros_manobra_adquire_gso(t, m, X)  # Atualiza variáveis globais

    # Derivada do vetor de estado
    Xp = [Vp, Ap, phip, rp, deltap, lonp]
    return Xp
