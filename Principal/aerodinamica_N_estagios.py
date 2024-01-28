from modelo_aerodinamico import modelo_aerodinamico
import parametros

def aerodinamica_N_estagios(t, V, h, M, Kn, T, rho, R):
    ts = parametros.ts
    Sr = parametros.Sr
    fator_de_correcao_arrasto = parametros.fator_correcao_arrasto
    # Modelo de arrasto conforme a referencia
    # TEWARI, A. Atmospheric and Space Flight Dynamics:
    # Modelling and simulation with MATLAB and Simulink. Boston: Birkhauser, 2007.
    # Exemplo 12.6
    # Vale para foguete de 1, 2 ou 3 estagios, com carga util
    # Assume-se o mesmo coeficiente de arrasto para o foguete todo e seus
    # estagios. A magnitude do arrasto eh alterada em funcao da area de
    # referencia de cada estagio

    # Coeficiente de arrasto em funcao do numero de Mach e de Knuden
    coeficiente_de_arrasto_machKnudesen = modelo_aerodinamico(V, h, M, Kn, T, R)
    # Fator de correcao do arrasto a partir de dados de tunel de vento
    coeficiente_de_arrasto_machKnudesen = fator_de_correcao_arrasto * coeficiente_de_arrasto_machKnudesen

    # A area de referencia depende do estagio atual
    # Numero de estagios
    N = len(ts)
    if N == 1:
        S = area1estagio(t, ts, Sr)
    elif N == 2:
        S = area2estagios(t, ts, Sr)
    elif N == 3:
        S = area3estagios(t, ts, Sr)

    # Forcas
    D = 0.5 * rho * V ** 2 * S * coeficiente_de_arrasto_machKnudesen
    fy = 0
    L = 0

    return D, fy, L


def area1estagio(t, ts, Sr):
    if t <= ts[0]:
        # Foguete e carga util
        S = Sr[0]
    else:
        S = Sr[1]  # Carga util
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