import numpy as np
from scipy.optimize import fsolve
from src.domain.utilidades_mecanica_orbital.Orbitas.ModeloOrbita import Orbita

def resolveEqKepler(t, orbita: Orbita):
    """
    Resolve a equação de Kepler para uma órbita elíptica, calculando a anomalia excêntrica para um dado tempo.

    Parâmetros:
    t (float): Tempo em segundos.
    orbita (Orbita): Objeto da classe Orbita contendo os elementos orbitais clássicos.

    Retorna:
    float: Anomalia excêntrica no instante t.
    """
    # Extrai os elementos orbitais do objeto Orbita
    a = orbita.semi_eixo_maior
    e = orbita.excentricidade
    mu = orbita.mu
    tau = orbita.tempo_de_periastro

    # Calcula o movimento médio
    n = np.sqrt(mu / a ** 3)

    # Calcula a anomalia média
    M = n * (t - tau)

    # Passa os parâmetros para a função objetivo
    dados = (e, M)

    # Chute inicial para a solução
    E0 = M

    # Resolve a equação de Kepler usando fsolve
    E = fsolve(Kepler, E0, args=dados)
    return E[0]

def propagaEliptica(t, orbita):
    """
    Propaga a órbita elíptica no referencial inercial para um tempo dado.

    Parâmetros:
    t (float): Tempo em segundos.
    orbita (Orbita): Objeto da classe Orbita contendo os elementos orbitais clássicos.

    Retorna:
    tuple: Anomalia verdadeira (theta), vetor posição (Ri), vetor velocidade (Vi).
    """
    # Extrai os elementos orbitais do objeto Orbita
    a = orbita.semi_eixo_maior
    e = orbita.excentricidade
    i = orbita.inclinacao
    omega = orbita.argumento_de_periastro
    OMEGA = orbita.longitude_celeste
    mu = orbita.mu
    tau = orbita.tempo_de_periastro

    # Resolve a equação de Kepler para obter a anomalia excêntrica
    E = resolveEqKepler(t, orbita)

    # Calcula a anomalia verdadeira
    theta = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E / 2), np.sqrt(1 - e) * np.cos(E / 2))

    # Calcula a distância radial
    r = a * (1 - e ** 2) / (1 + e * np.cos(theta))

    # Vetor posição no referencial perifocal
    R_perifocal = np.array([r * np.cos(theta), r * np.sin(theta), 0])

    # Vetor velocidade no referencial perifocal
    V_perifocal = np.array([-np.sin(E), np.sqrt(1 - e ** 2) * np.cos(E), 0]) * np.sqrt(mu * a) / r

    # Converte para o referencial inercial
    R_inercial = perifocal_to_inercial(R_perifocal, OMEGA, i, omega)
    V_inercial = perifocal_to_inercial(V_perifocal, OMEGA, i, omega)

    return theta, R_inercial, V_inercial

def perifocal_to_inercial(vector, OMEGA, i, omega):
    """
    Converte um vetor do referencial perifocal para o referencial inercial.

    Parâmetros:
    vector (np.array): Vetor no referencial perifocal.
    OMEGA (float): Longitude do nodo ascendente em radianos.
    i (float): Inclinação da órbita em radianos.
    omega (float): Argumento do periastro em radianos.

    Retorna:
    np.array: Vetor no referencial inercial.
    """
    # Matriz de rotação composta
    rotation_matrix = np.array([
        [np.cos(OMEGA) * np.cos(omega) - np.sin(OMEGA) * np.sin(omega) * np.cos(i),
         -np.cos(OMEGA) * np.sin(omega) - np.sin(OMEGA) * np.cos(omega) * np.cos(i),
         np.sin(OMEGA) * np.sin(i)],
        [np.sin(OMEGA) * np.cos(omega) + np.cos(OMEGA) * np.sin(omega) * np.cos(i),
         -np.sin(OMEGA) * np.sin(omega) + np.cos(OMEGA) * np.cos(omega) * np.cos(i),
         -np.cos(OMEGA) * np.sin(i)],
        [np.sin(omega) * np.sin(i),
         np.cos(omega) * np.sin(i),
         np.cos(i)]
    ])

    # Aplica a matriz de rotação ao vetor
    return np.dot(rotation_matrix, vector)

def Kepler(E, dados):
    """
    Função objetivo para resolver a equação de Kepler.

    Parâmetros:
    E (float): Anomalia excêntrica.
    dados (tuple): Contém a excentricidade (e) e a anomalia média (M).

    Retorna:
    float: Diferença entre os lados da equação de Kepler.
    """
    e, M = dados
    return E - e * np.sin(E) - M

def KeplerHiperbolica(H, dados):
    """
    Função objetivo para resolver a equação de Kepler para órbitas hiperbólicas.

    Parâmetros:
    H (float): Anomalia hiperbólica.
    dados (tuple): Contém a excentricidade (e) e a anomalia média hiperbólica (Mh).

    Retorna:
    float: Diferença entre os lados da equação de Kepler hiperbólica.
    """
    e, Mh = dados
    return e * np.sinh(H) - H - Mh

def resolveEqKeplerHiperbolica(tempo: float, orbita: Orbita):
    """
    Resolve a equação de Kepler para uma órbita hiperbólica, calculando a anomalia hiperbólica para um dado tempo.

    Parâmetros:
    tempo (float): Tempo em segundos.
    orbita (Orbita): Objeto da classe Orbita contendo os elementos orbitais clássicos.

    Retorna:
    float: Anomalia hiperbólica no instante tempo.
    """
    tau = orbita.tempo_de_periastro
    p = orbita.calcular_parametro_orbital()
    e = orbita.excentricidade
    Mh = (e ** 2 - 1) ** (3 / 2) * np.sqrt(orbita.mu / p ** 3) * (tempo - tau)
    dados = (e, Mh)
    H0 = Mh
    H = fsolve(KeplerHiperbolica, H0, args=dados)
    return H

def matrizTransicaoEstado(theta, theta0, orbita):
    """
    Determina a matriz de transição de estado cujos elementos são os coeficientes de Lagrange.

    Parâmetros:
    theta (float): Anomalia verdadeira em radianos.
    theta0 (float): Anomalia verdadeira inicial em radianos.
    orbita (Orbita): Objeto da classe Orbita contendo os elementos orbitais clássicos.

    Retorna:
    np.array: Matriz de transição de estado (2x2).
    """
    e = orbita.excentricidade
    mu = orbita.mu
    p = orbita.calcular_parametro_orbital()
    h = np.sqrt(p * mu)
    r0 = p / (1 + e * np.cos(theta0))
    r = p / (1 + e * np.cos(theta))
    f = 1 + (r / p) * (np.cos(theta - theta0) - 1)
    g = (r * r0 / h) * np.sin(theta - theta0)
    fp = -(h / p ** 2) * (np.sin(theta - theta0) + e * (np.sin(theta) - np.sin(theta0)))
    gp = 1 + (r0 / p) * (np.cos(theta - theta0) - 1)
    return np.array([[f, g], [fp, gp]])
