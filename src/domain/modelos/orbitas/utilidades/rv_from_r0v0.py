import numpy as np

from src.domain.modelos.orbitas.utilidades.coeficientes_de_lagrange import f_and_g, fDot_and_gDot
from src.domain.modelos.orbitas.utilidades.kepler_universal import kepler_U


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

    # Magnitudes of R0 and V0
    r0 = np.linalg.norm(R0)
    v0 = np.linalg.norm(V0)

    # Initial radial velocity
    vr0 = np.dot(R0, V0) / r0

    # Reciprocal of the semimajor axis (from the energy equation)
    alpha = 2 / r0 - v0 ** 2 / mu

    # Compute the universal anomaly
    x = kepler_U(delta_t, r0, vr0, alpha, mu)

    # Compute the f and g functions
    f, g = f_and_g(x, delta_t, r0, alpha, mu)

    # Compute the final position vector
    R = f * R0 + g * V0

    # Compute the magnitude of R
    r = np.linalg.norm(R)

    # Compute the derivatives of f and g
    fdot, gdot = fDot_and_gDot(x, r, r0, alpha, mu)

    # Compute the final velocity
    V = fdot * R0 + gdot * V0

    return R, V








