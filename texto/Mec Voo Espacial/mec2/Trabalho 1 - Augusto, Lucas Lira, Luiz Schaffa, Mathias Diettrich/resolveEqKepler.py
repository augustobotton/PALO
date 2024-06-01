import parametros
import numpy as np
from scipy.optimize import fsolve
#
# Órbita eliptica: calculo da propagacao temporal da anomalia excentrica
#
def Kepler(E,*dados):
# Funcao objetivo: Equacao de Kepler
# Entradas
# E: [rad] anomalia excentrica - incógnita
# *dados: parâmetros da função
# Saida
# y: quando igual a zero, a equacao de Kepler esta resolvida
#
# Recebimento de parametros
    e,M=dados
# Quando y=0, a equacao de Kepler esta resolvida
    y=E-e*np.sin(E)-M
    return y
#
def resolveEqKepler(t):
# Órbita eliptica: calcula a anomalia excentrica para dado tempo
#
# Entrada
# t: [s] tempo
# Saida
# E: [rad] anomalia excentrica
# Recebe dados pelo arquivo parametros
    tau = parametros.TAU # tempo de periastro
    a = parametros.A # Semi eixo maior
    mu = parametros.MU # Constante gravitacional do planeta
    e = parametros.E # Excentricidade da orbita
# Movimento medio
    n = np.sqrt(mu/a**3)
# Anomalia media
    M = n*(t-tau)
# Passagem de parâmetros para a função objetivo
    dados = (e,M)
# Chute inicial
    E0 = M # Seria o resultado da orbita circular
# Resolve numericamente a equacao de Kepler com a funcao fsolve do pacote scipy.optimize
    E = fsolve(Kepler,E0,args=dados)
    return E
