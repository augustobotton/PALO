import numpy as np

from src.domain.utilidades_mecanica_orbital.Utilidades.CoeficientesDeLagrange import f_and_g, fDot_and_gDot
from src.domain.utilidades_mecanica_orbital.Utilidades.KeplerUniversal import kepler_U


def rv_from_r0v0(R0, V0, delta_t, mu):
    """
    Calcula o vetor de posição final e o vetor de velocidade final dado um vetor de posição inicial,
    um vetor de velocidade inicial e o tempo decorrido. Baseado na solução da equação de Kepler Universal.

    Parâmetros:
        R0 (np.array): Vetor de posição inicial (km)
        V0 (np.array): Vetor de velocidade inicial (km/s)
        t (float): Tempo decorrido (s)

    Retorna:
        np.array: Vetor de posição final (km)
        np.array: Vetor de velocidade final (km/s)
    """

    r0 = np.linalg.norm(R0)
    v0 = np.linalg.norm(V0)
    vr0 = np.dot(R0, V0) / r0
    alpha = 2 / r0 - v0 ** 2 / mu

    x = kepler_U(delta_t, r0, vr0, alpha, mu)
    f, g = f_and_g(x, delta_t, r0, alpha, mu)
    R = f * R0 + g * V0
    r = np.linalg.norm(R)
    fdot, gdot = fDot_and_gDot(x, r, r0, alpha,mu)
    V = fdot * R0 + gdot * V0

    return R, V









