import modelo_aerodinamico
from dinamica_foguete import *


def area1estagio(t, ts, Sr):
    if t <= ts[0]:
        # Foguete e carga útil
        S = Sr[0]
    else:
        S = Sr[1]  # Carga útil
    return S

def area2estagios(t, ts, Sr):
    if t <= ts[0]:
        # Todos os estágios
        S = Sr[0]
    elif t <= ts[1]:
        # Segundo estágio e carga útil
        S = Sr[1]
    else:
        # Carga útil
        S = Sr[2]
    return S

def area3estagios(t, ts, Sr):
    if t <= ts[0]:
        # Todos os estágios
        S = Sr[0]
    elif t <= ts[1]:
        # Segundo estágio
        S = Sr[1]
    elif t <= ts[2]:
        # Terceiro estágio e carga útil
        S = Sr[2]
    else:
        # Carga útil
        S = Sr[3]
    return S

def aerodinamica_N_estagios(t, V, h, M, Kn, T, R, ts, Sr, fc):
    global rho
    # Coeficiente de arrasto em função do número de Mach e de Knuden
    CD = modelo_aerodinamico(V, h, M, Kn, T, R)
    # Fator de correção do arrasto a partir de dados de túnel de vento
    CD = fc * CD
    
    # A área de referência depende do estágio atual
    # Número de estágios
    N = len(ts)
    if N == 1:
        S = area1estagio(t, ts, Sr)
    elif N == 2:
        S = area2estagios(t, ts, Sr)
    elif N == 3:
        S = area3estagios(t, ts, Sr)
    
    # Forças
    D = 0.5 * rho * V**2 * S * CD
    fy, L = 0, 0
    
    return D, fy, L