import numpy as np
from scipy.optimize import fsolve

def calcular_massas_estagios(isp, forca_propulsiva_kN, razao_estrutural, V_bo, chute_inicial_eta, m_PL):
    def equacao_para_resolver(eta, c, razao_estrutural, V_bo):
        termo1 = np.sum(c * np.log(c * eta - 1))
        termo2 = np.log(eta) * np.sum(c)
        termo3 = np.sum(c * np.log(c * razao_estrutural))
        return termo1 - termo2 - termo3 - V_bo

    def resolver_para_eta(c, razao_estrutural, V_bo, chute_inicial_eta):
        solucao_eta = fsolve(equacao_para_resolver, chute_inicial_eta, args=(c, razao_estrutural, V_bo))
        return solucao_eta[0]

    def calcula_duracao_queima(massa_propelente, isp, F):
        return massa_propelente * isp / F

    # Convertendo Isp para c (km/s)
    g0 = 9.80665  # m/s^2
    c = isp * g0 / 1000  # Convertendo para km/s

    N = len(c)  # Número de estágios baseado no tamanho do vetor c
    solucao_eta = resolver_para_eta(c, razao_estrutural, V_bo, chute_inicial_eta)
    print("Solução para eta:", solucao_eta)
    n = (c * solucao_eta - 1) / (c * razao_estrutural * solucao_eta)

    m = np.zeros(N)  # Massas dos estágios
    mE = np.zeros(N)  # Massas vazias dos estágios
    mp = np.zeros(N)  # Massas de propelente dos estágios

    # Calcular a massa da última etapa (topo)
    m[N-1] = (n[N-1] - 1) / (1 - n[N-1] * razao_estrutural[N-1]) * m_PL

    # Calcular as massas para as demais etapas de maneira reversa
    for i in range(N-2, -1, -1):
        m[i] = (n[i] - 1) / (1 - n[i] * razao_estrutural[i]) * (np.sum(m[i+1:]) + m_PL)

    # Calcular as massas vazias e de propelente para cada estágio
    for i in range(N):
        mE[i] = razao_estrutural[i] * m[i]
        mp[i] = m[i] - mE[i]

    # Converter força propulsiva de kN para N
    forca_propulsiva_N = forca_propulsiva_kN * 1000

    # Calcular a duração da queima para cada estágio
    duracao_queima = calcula_duracao_queima(mp, isp, forca_propulsiva_N)

    return m, mE, mp, duracao_queima

# Exemplo de uso
isp = np.array([290, 420, 420])  # Impulso específico em segundos
forca_propulsiva_kN = np.array([20040, 4450, 890])  # Força propulsiva em kN
razao_estrutural = np.array([0.1, 0.15, 0.2])  # Valores exemplo para epsilon_i
V_bo = 9.5
chute_inicial_eta = 0.5
m_PL = 85000  # Massa da carga útil

massas_estagios, massas_vazias, massas_propelente, duracao_queima = calcular_massas_estagios(isp, forca_propulsiva_kN, razao_estrutural, V_bo, chute_inicial_eta, m_PL)
print("Massas dos estágios:", massas_estagios)
print("Massas vazias dos estágios:", massas_vazias)
print("Massas de propelente dos estágios:", massas_propelente)
print("Duração de queima dos estágios (em segundos):", duracao_queima)
