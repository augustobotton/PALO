import numpy as np

from src.domain.modelos.orbitas.utilidades.rv_from_r0v0 import rv_from_r0v0


def dinamica_perturbada_J2(t, f, *args):
    """
    Calcula a dinâmica de uma órbita perturbada pelo efeito J2 da Terra.

    :param t: Tempo atual da integração (em segundos).
    :param f: Vetor de estado atual (6x1) contendo posição e velocidade perturbadas.
    :param args: Argumentos adicionais:
        - mu: Parâmetro gravitacional da Terra (km^3/s^2).
        - R0: Vetor de posição inicial (km).
        - V0: Vetor de velocidade inicial (km/s).
        - t0: Tempo inicial da integração (s).
        - RE: Raio equatorial da Terra (km).
        - J2: Coeficiente J2 da Terra.
    :return: Vetor concatenado de derivadas de posição e velocidade (6x1).
    """
    mu = args[0]
    R0 = args[1]
    V0 = args[2]
    t0 = args[3]
    RE = args[4]
    J2 = args[5]

    # Desempacota o vetor de estado
    del_r = f[:3]
    del_v = f[3:6]

    # Computa o vetor de estado na órbita osculadora no tempo t
    Rosc, Vosc = rv_from_r0v0(R0, V0, t - t0, mu)

    # Calcula os componentes do vetor de estado na órbita perturbada
    Rpp = Rosc + del_r
    Vpp = Vosc + del_v
    rosc = np.linalg.norm(Rosc)
    rpp = np.linalg.norm(Rpp)

    # Computa a aceleração perturbativa devido ao efeito J2
    xx = Rpp[0]
    yy = Rpp[1]
    zz = Rpp[2]
    fac = 1.5 * J2 * (mu / rpp ** 2) * (RE / rpp) ** 2
    ap = -fac * np.array([
        (1 - 5 * (zz / rpp) ** 2) * (xx / rpp),
        (1 - 5 * (zz / rpp) ** 2) * (yy / rpp),
        (3 - 5 * (zz / rpp) ** 2) * (zz / rpp)
    ])

    # Computa a aceleração perturbativa total
    F = 1 - (rosc / rpp) ** 3
    del_a = (-mu / rosc ** 3) * (del_r - F * Rpp) + ap

    return np.concatenate([del_v, del_a])
