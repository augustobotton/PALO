import numpy as np

from src.domain.utilidades_mecanica_orbital.Utilidades.funcoes_stump import stumpC, stumpS


def f_and_g(x, t, ro, a, mu):
    """
    Calcula os coeficientes de Lagrange f e g dados os parâmetros de entrada.

    Parâmetros:
        x (float): Anomalia universal após o tempo t (km^0.5)
        t (float): Tempo decorrido desde ro (s)
        ro (float): Posição radial no tempo t0 (km)
        a (float): Recíproco do semieixo maior (1/km)

    Retorna:
        float: Coeficiente de Lagrange f (adimensional)
        float: Coeficiente de Lagrange g (s)
    """
    z = a * x ** 2
    C = stumpC(z)
    S = stumpS(z)
    f = 1 - x ** 2 / ro * C
    g = t - 1 / np.sqrt(mu) * x ** 3 * S

    return f, g


def fDot_and_gDot(x, r, ro, a, mu):
    """
    Calcula as derivadas temporais dos coeficientes de Lagrange f e g.

    Parâmetros:
        x (float): Anomalia universal após o tempo t (km^0.5)
        r (float): Posição radial após o tempo t (km)
        ro (float): Posição radial no tempo t0 (km)
        a (float): Recíproco do semieixo maior (1/km)

    Retorna:
        float: Derivada temporal do coeficiente de Lagrange f (1/s)
        float: Derivada temporal do coeficiente de Lagrange g (adimensional)
    """
    z = a * x ** 2
    S = stumpS(z)
    C = stumpC(z)
    fdot = np.sqrt(mu) / (r * ro) * (z * S - 1) * x
    gdot = 1 - x ** 2 / r * C

    return fdot, gdot
