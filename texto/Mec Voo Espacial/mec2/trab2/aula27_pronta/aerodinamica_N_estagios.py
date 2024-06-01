import numpy as np
from modelo_aerodinamico import modelo_aerodinamico


def aerodinamica_N_estagios(t, V, h, M, Kn, T, rho, R,fc,ts,Sr):
    
    # Area de referencia de foguete de 1 estagio
    def area1estagio(t, ts, Sr):
        if t <= ts[0]:
            # Foguete e carga util
            S = Sr[0]
        else:
            S = Sr[1]   # Carga util
        return S

# Area de referencia de foguete de 2 estagios
    def area2estagios(t, ts, Sr):
        if t <= ts[0]:
            # Todos os estagios
            S = Sr[0]
        elif t <= ts[1]:
            # Segundo estagio e carga util
            S = Sr[1]
        else:
            # Carga util
            S = Sr[2]
        return S

    # Area de referencia de foguete de 3 estagios
    def area3estagios(t, ts, Sr):
        if t <= ts[0]:
            # Todos os estagios
            S = Sr[0]
        elif t <= ts[1]:
            # Segundo estagio
            S = Sr[1]
        elif t <= ts[2]:
            # Terceiro estagio e carga util
            S = Sr[2]
        else:
            # Carga util
            S = Sr[3]
        return S



    # Coeficiente de arrasto em funcao do numero de Mach e de Knuden
    CD = modelo_aerodinamico(V, h, M, Kn, T, R)
    # Fator de correcao do arrasto a partir de dados de tunel de vento
    CD = fc * CD
    # A area de referencia depende do estagio atual
    # Numero de estagios
    N = len(ts)
    if N == 1:
        S = area1estagio(t, ts, Sr)
    elif N == 2:
        S = area2estagios(t, ts, Sr)
    else: 
        S = area3estagios(t, ts, Sr)
        
    # Forcas
    D = 0.5 * rho * V**2 * S * CD
    fy = 0
    L = 0
    return D, fy, L

