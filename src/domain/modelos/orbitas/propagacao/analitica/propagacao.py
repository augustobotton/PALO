import numpy as np
from scipy.optimize import fsolve

from src.domain.modelos.orbitas.Orbita import Orbita
from src.domain.modelos.orbitas.utilidades.funcoes_conversao import perifocal_para_inercial


def resolveEqKepler(t, orbita: Orbita):
    """
    Resolve a equação de Kepler para uma órbita elíptica, calculando a anomalia excêntrica.

    Parâmetros:
    t (float): Tempo em segundos.
    orbita (Orbita): Objeto da classe Orbita com elementos orbitais clássicos.

    Retorna:
    float: Anomalia excêntrica no instante t.
    """
    a, e, mu, tau = orbita.semi_eixo_maior, orbita.excentricidade, orbita.mu, orbita.tempo_de_periastro
    n = np.sqrt(mu / a ** 3)
    M = n * (t - tau)
    E = fsolve(lambda E: E - e * np.sin(E) - M, M)
    return E[0]


def propagaEliptica(t, orbita):
    """
    Propaga a órbita elíptica no referencial inercial para um tempo dado.

    Parâmetros:
    t (float): Tempo em segundos.
    orbita (Orbita): Objeto da classe Orbita com elementos orbitais clássicos.

    Retorna:
    tuple: Anomalia verdadeira (theta), vetor posição (Ri), vetor velocidade (Vi).
    """
    E = resolveEqKepler(t, orbita)
    theta = 2 * np.arctan2(np.sqrt(1 + orbita.excentricidade) * np.sin(E / 2),
                           np.sqrt(1 - orbita.excentricidade) * np.cos(E / 2))
    r = orbita.semi_eixo_maior * (1 - orbita.excentricidade ** 2) / (1 + orbita.excentricidade * np.cos(theta))
    R_perifocal = np.array([r * np.cos(theta), r * np.sin(theta), 0])
    V_perifocal = np.array([-np.sin(E), np.sqrt(1 - orbita.excentricidade ** 2) * np.cos(E), 0]) * np.sqrt(
        orbita.mu * orbita.semi_eixo_maior) / r
    R_inercial, V_inercial = perifocal_para_inercial(R_perifocal, orbita.longitude_celeste, orbita.inclinacao,
                                                     orbita.argumento_de_periastro), perifocal_para_inercial(
        V_perifocal, orbita.longitude_celeste, orbita.inclinacao, orbita.argumento_de_periastro)
    return theta, R_inercial, V_inercial


def resolveEqKeplerHiperbolica(t, orbita: Orbita):
    """
    Resolve a equação de Kepler para uma órbita hiperbólica, calculando a anomalia hiperbólica.

    Parâmetros:
    t (float): Tempo em segundos.
    orbita (Orbita): Objeto da classe Orbita com elementos orbitais clássicos.

    Retorna:
    float: Anomalia hiperbólica no instante t.
    """
    tau, p, e = orbita.tempo_de_periastro, orbita.calcular_parametro_orbital(), orbita.excentricidade
    Mh = (e ** 2 - 1) ** 1.5 * np.sqrt(orbita.mu / p ** 3) * (t - tau)
    H = fsolve(lambda H: e * np.sinh(H) - H - Mh, Mh)
    return H[0]


def matrizTransicaoEstado(theta, theta0, orbita):
    """
    Determina a matriz de transição de estado com coeficientes de Lagrange.

    Parâmetros:
    theta (float): Anomalia verdadeira em radianos.
    theta0 (float): Anomalia verdadeira inicial em radianos.
    orbita (Orbita): Objeto da classe Orbita com elementos orbitais clássicos.

    Retorna:
    np.array: Matriz de transição de estado (2x2).
    """
    e, mu, p = orbita.excentricidade, orbita.mu, orbita.calcular_parametro_orbital()
    h = np.sqrt(p * mu)
    r0 = p / (1 + e * np.cos(theta0))
    r = p / (1 + e * np.cos(theta))
    f = 1 + (r / p) * (np.cos(theta - theta0) - 1)
    g = r * r0 / h * np.sin(theta - theta0)
    fp = -(h / p ** 2) * (np.sin(theta - theta0) + e * (np.sin(theta) - np.sin(theta0)))
    gp = 1 + (r0 / p) * (np.cos(theta - theta0) - 1)
    return np.array([[f, g], [fp, gp]])
