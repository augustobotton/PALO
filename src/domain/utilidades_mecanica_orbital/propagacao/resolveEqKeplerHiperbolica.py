import numpy as np
from scipy.optimize import fsolve
#
# Órbita hiperbólica: cálculo da propagacao temporal da anomalia hiperbólica
#
def KeplerHiperbolica(H,*dados):
    # Função objetivo: equação de Kepler hiperbólica
    # Entrada
    # H: [rad] anomalia hiperbólica
    # *dados: parâmetros da função
    # Saida
    # y: quando igual a zero, a equacao de Kepler hiperbólica está resolvida
    #
    # Recebimento de parametros
    e,Mh=dados
    # Quando y=0, a equacao de Kepler hiperbólica está resolvida
    y = e*np.sinh(H)-H-Mh
    return y
    #

def resolveEqKeplerHiperbolica(t):
    # Órbita hiperbólica: calcula a anomalia hiperbólica para dado tempo
    # Entrada
    # t: [s] tempo
    # Saida
    #
    # H: [rad] anomalia hiperbólica
    # Recebe dados pelo arquivo parametros
    tau = parametros.TAU # tempo de periastro
    p = parametros.P
    # Parâmetro
    mu = parametros.MU # Constante gravitacional do planeta
    e = parametros.E # Excentricidade da orbita
    # Anomalia media hiperbolica
    Mh=(e**2-1)**(3/2)*np.sqrt(mu/p**3)*(t-tau)
    # Passagem de parâmetros para a função objetivo
    dados = (e,Mh)
    # Chute inicial
    H0 = Mh

    # Resolve numericamente a equacao de Kepler hiperbólica com a funcao fsolve do
    # pacote scipy.optimize
    H = fsolve(KeplerHiperbolica,H0,args=dados)
    return H
