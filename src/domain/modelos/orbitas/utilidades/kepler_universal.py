import numpy as np

from src.domain.modelos.orbitas.utilidades.funcoes_stump import stumpC, stumpS


def kepler_U(dt, ro, vro, a, mu):
    """
    Resolve a equação de Kepler universal para a anomalia universal usando o método de Newton.

    Parâmetros:
        dt (float): Tempo desde x = 0 (s)
        ro (float): Posição radial (km) quando x = 0
        vro (float): Velocidade radial (km/s) quando x = 0
        a (float): Recíproco do semieixo maior (1/km)

    Retorna:
        float: A anomalia universal (km^0.5)
    """
    # Set an error tolerance and a limit on the number of iterations
    error = 1.e-8
    nMax = 1000

    # Starting value for x
    x = np.sqrt(mu) * np.abs(a) * dt

    # Iterate on Equation 3.65 until convergence occurs within the error tolerance
    n = 0
    ratio = 1

    while np.abs(ratio) > error and n <= nMax:
        n += 1
        C = stumpC(a * x ** 2)
        S = stumpS(a * x ** 2)
        F = ro * vro / np.sqrt(mu) * x ** 2 * C + (1 - a * ro) * x ** 3 * S + ro * x - np.sqrt(mu) * dt
        dFdx = ro * vro / np.sqrt(mu) * x * (1 - a * x ** 2 * S) + (1 - a * ro) * x ** 2 * C + ro
        ratio = F / dFdx
        x -= ratio

    # Deliver a value for x, but report that nMax was reached
    if n > nMax:
        print(f"\n **No. iterations of Kepler’s equation = {n}")
        print(f"\n F/dFdx = {ratio} \n")

    return x
