from grav_axisimetrico import *
import numpy as np
from scipy.constants import g
from parametros_manobra_adquire_gso import parametros_manobra_adquire_gso
from aerodinamica_N_estagios import aerodinamica_N_estagios
from parametrosglobais import *
from propulsao_N_estagios import propulsao_N_estagios
from atm_padrao import atm_padrao
# As funções "propulsao_N_estagios", "atm_padrao", "aerodinamica_N_estagios", "grav_axisimetrico" 
# e "parametros_manobra_adquire_gso" precisam ser definidas. Os nomes dessas funções são traduções diretas do MATLAB, mas 
# as implementações dessas funções dependem do que elas estão fazendo exatamente no seu código original.

def dinamica_foguete(t, X):
    global we, Requat, lc, dT, h0, l_trilho, rho
    # Vetor de estado
    V, A, phi, r, delta = X[0], X[1], X[2], X[3], X[4]
    
    if V < 0:
        V = 0  # Evita velocidade negativa 

    # Função para cálculo da massa e tração em função do tempo
    ft, m, mu, epsl = propulsao_N_estagios(t, X)
    
    # Função para cálculo do modelo atmosférico
    h = r - Requat  # Altitude
    T, _, rho, M, _, Kn, _, _, R = atm_padrao(h, V, lc, dT)
    
    # Cálculo do modelo aerodinâmico
    D, fy, L = aerodinamica_N_estagios(t, V, h, M, Kn, T, rho, R)
    
    # Cálculo da gravidade
    gc, gd = grav_axisimetrico(r, delta)

    # Equações de cinemática de translação
    rp = V * np.sin(phi)
    deltap = (V / r) * np.cos(phi) * np.cos(A)
    lonp = (V * np.cos(phi) * np.sin(A)) / (r * np.cos(delta))

    # Equações de dinâmica de translação
    Vp = (1 / m) * (ft * np.cos(epsl) * np.cos(mu) - D - m * gc * np.sin(phi) + m * gd * np.cos(phi) * np.cos(A) - 
        m * we**2 * r * np.cos(delta) * (np.cos(phi) * np.cos(A) * np.sin(delta) - np.sin(phi) * np.cos(delta)))
    Ap = (1 / (m * V * np.cos(phi))) * (m * (V**2 / r) * np.cos(phi)**2 * np.sin(A) * np.tan(delta) + ft * np.sin(mu) + fy 
        - m * gd * np.sin(A) + m * we**2 * r * np.sin(A) * np.sin(delta) * np.cos(delta) - 2 * m * we * V * 
        (np.sin(phi) * np.cos(A) * np.cos(delta) - np.cos(phi) * np.sin(delta)))
    phip = (1 / (m * V)) * (m * (V**2 / r) * np.cos(phi) + ft * np.sin(epsl) * np.cos(mu) + L - m * gc * np.cos(phi) - 
        m * gd * np.sin(phi) * np.cos(A) + m * we**2 * r * np.cos(delta) * (np.sin(phi) * np.cos(A) * np.sin(delta) + 
        np.cos(phi) * np.cos(delta)) + 2 * m * we * V * np.sin(A) * np.cos(delta))

    # Saturação da altitude
    if h < 0:  # Altitude negativa não é permitida
        # Mantém as derivadas nulas
        rp = 0; deltap = 0; lonp = 0; Vp = 0; Ap = 0; phip = 0
    
    # Modela o trilho de lançamento
    H = h - h0  # Altura
    if (H <= l_trilho) and (t <= 10):  # Verifica se a altura é menor que l_trilho nos primeiros segundos da simulação
        Ap = 0; phip = 0  # Anula as derivadas dos ângulos de orientação da velocidade

    # Modela a navegação para o segundo disparo do motor do terceiro estágio
    parametros_manobra_adquire_gso(t, m, X)  # Atualiza variáveis globais

    # Derivada do vetor de estado
    Xp = np.array([Vp, Ap, phip, rp, deltap, lonp])
  
    return Xp
