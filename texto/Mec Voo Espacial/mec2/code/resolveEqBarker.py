import parametros
import numpy as np
#
def resolveEqBarker(t):
# Orbita parabolica: Calcula a anomalia verdadeira dado o tempo
# Entrada
# t: [s] tempo
# Saida
# theta: [rad] anomalia verdadeira
#
# Recebe dados pelo arquivo parametros
    tau = parametros.TAU # tempo de periastro
    p = parametros.P
# Parametro
    mu = parametros.MU # Constante gravitacional do planeta
# Anomalia media parabolica
    Mp = np.sqrt(mu/p**3)*(t-tau)
# Solucao analitica - tangente de meio theta
    tantheta2 = (3*Mp+np.sqrt(1+9*Mp**2))**(1/3)-(3*Mp+np.sqrt(1+9*Mp**2))**(-1/3)
# Anomalia verdadeira
    theta = 2*np.arctan(tantheta2)
    return theta
