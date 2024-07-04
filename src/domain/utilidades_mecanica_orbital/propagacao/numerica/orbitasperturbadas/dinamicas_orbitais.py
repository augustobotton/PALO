import numpy as np

from src.domain.utilidades_mecanica_orbital.orbitalUtils.rv_from_r0v0 import rv_from_r0v0


def aceleracao_perturbativa_terceiro_corpo(
        posicao_m2: np.ndarray,
        posicao_m3: np.ndarray,
        parametro_gravitacional_m3: float
) -> np.ndarray:
    """
    Calcula a aceleração perturbativa de um terceiro corpo sobre uma massa de prova em órbita.

    :param posicao_m2: Vetor 3x1 da posição da massa de prova m2 com respeito ao primário m1 (m).
    :param posicao_m3: Vetor 3x1 da posição do terceiro corpo m3 com respeito ao primário m1 (m).
    :param parametro_gravitacional_m3: Parâmetro gravitacional do terceiro corpo m3 (m^3/s^2).
    :return: Vetor 3x1 da aceleração perturbativa do terceiro corpo sobre m2 (m/s^2).
    """
    # Posição de m2 com respeito a m3
    posicao_relativa_m2_m3 = posicao_m2 - posicao_m3

    # Distâncias
    distancia_m2_m3 = np.linalg.norm(posicao_relativa_m2_m3)
    distancia_m1_m3 = np.linalg.norm(posicao_m3)

    # Aceleração perturbativa
    aceleracao_perturbativa = -parametro_gravitacional_m3 * (
            posicao_relativa_m2_m3 / distancia_m2_m3 ** 3 + posicao_m3 / distancia_m1_m3 ** 3
    )

    return aceleracao_perturbativa


def dinPert3corpo(t, X, orbita):
    """
    Entradas:
    t: [s] - tempo
    X: vetor de estado, 6x1
    Saída
    xp: derivada de x
    """
    # Posição e velocidade
    R = np.array(X[0:3])
    V = np.array(X[3:6])
    #Distância
    r = np.linalg.norm(R)
    # Posição do terceiro corpo com respeito ao primário
    _, R31, _ = propagacao.propagaEliptica(t, mu, orbita)
    # Aceleração perturbativa
    ad = aceleracao_perturbativa_terceiro_corpo(R, R31, mu3)
    # Cinemática
    Rp = V
    # Dinâmica
    Vp = -mu * R / r ** 3 + ad
    #Saída
    Xp = np.concatenate((Rp, Vp))
    return Xp


def dinamica_J2(t, f, *args):
    mu = args[0]
    R0 = args[1]
    V0 = args[2]
    t0 = args[3]
    RE = args[4]
    J2 = args[5]

    # Unpack the state vector
    del_r = f[:3]
    del_v = f[3:6]
    # Compute the state vector on the osculating orbit at time t
    Rosc, Vosc = rv_from_r0v0(R0, V0, t - t0, mu)
    # Calculate the components of the state vector on the perturbed orbit
    Rpp = Rosc + del_r
    Vpp = Vosc + del_v
    rosc = np.linalg.norm(Rosc)
    rpp = np.linalg.norm(Rpp)

    # Compute the J2 perturbing acceleration
    xx=Rpp[0]
    yy=Rpp[1]
    zz =Rpp[2]
    fac = 1.5 * J2 * (mu / rpp ** 2) * (RE / rpp) ** 2
    ap = -fac * np.array([
        (1 - 5 * (zz / rpp) ** 2) * (xx / rpp),
        (1 - 5 * (zz / rpp) ** 2) * (yy / rpp),
        (3 - 5 * (zz / rpp) ** 2) * (zz / rpp)
    ])

    # Compute the total perturbing acceleration
    F = 1 - (rosc / rpp) ** 3
    del_a = (-mu / rosc ** 3) * (del_r - F * Rpp) + ap

    return np.concatenate([del_v, del_a])


def dinamica_lambert(t, f, mu):
    x, y, z, vx, vy, vz = f
    r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
    ax = -mu * x / r ** 3
    ay = -mu * y / r ** 3
    az = -mu * z / r ** 3
    return [vx, vy, vz, ax, ay, az]
