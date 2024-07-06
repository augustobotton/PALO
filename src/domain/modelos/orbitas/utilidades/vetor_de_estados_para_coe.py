import numpy as np
# Usa Algoritimo 4.1 de Orbital Mechanics for Engineering Students  (Curtis, 2020)
def vetor_de_estados_para_coe(R, V, mu):
    """
    Calcula os elementos orbitais clássicos a partir do vetor de estado (R, V).

    Args:
        R (np.ndarray): Vetor de posição no espaço tridimensional.
        V (np.ndarray): Vetor de velocidade no espaço tridimensional.
        mu (float): Parâmetro gravitacional padrão do corpo central.

    Returns:
        tuple: Uma tupla contendo:
            - a (float): Semieixo maior da órbita.
            - e (float): Excentricidade da órbita.
            - incl (float): Inclinação da órbita em radianos.
            - RA (float): Ascensão reta do nó ascendente em radianos.
            - w (float): Argumento do periastro em radianos.
            - TA (float): Anomalia verdadeira em radianos.
            - h (float): Momento angular específico.
    """
    eps = 1.e-10

    # Calcula as normas dos vetores de posição e velocidade
    r = np.linalg.norm(R)
    v = np.linalg.norm(V)

    # Velocidade radial
    vr = np.dot(R, V) / r

    # Momento angular específico
    H = np.cross(R, V)
    h = np.linalg.norm(H)

    # Inclinação
    incl = np.arccos(H[2] / h)

    # Nó ascendente
    K = np.array([0, 0, 1])
    N = np.cross(K, H)
    n = np.linalg.norm(N)

    # Ascensão reta do nó ascendente
    if n != 0:
        RA = np.arccos(N[0] / n)
        if N[1] < 0:
            RA = 2 * np.pi - RA
    else:
        RA = 0

    # Vetor excentricidade
    E = (1 / mu) * ((v ** 2 - (mu / r)) * R - (r * vr * V))
    e = np.linalg.norm(E)

    # Argumento do periastro
    if n != 0:
        if e > eps:
            w = np.arccos(np.dot(N, E) / (n * e))
            if E[2] < 0:
                w = 2 * np.pi - w
        else:
            w = 0
    else:
        w = 0

    # Anomalia verdadeira
    if e > eps:
        TA = np.arccos(np.dot(E, R) / (e * r))
        if vr < 0:
            TA = 2 * np.pi - TA
    else:
        cp = np.cross(N, R)
        if cp[2] >= 0:
            TA = np.arccos(np.dot(N, R) / (n * r))
        else:
            TA = 2 * np.pi - np.arccos(np.dot(N, R) / (n * r))

    # Semieixo maior
    a = h ** 2 / mu / (1 - e ** 2)

    return a, e, incl, RA, w, TA, h
