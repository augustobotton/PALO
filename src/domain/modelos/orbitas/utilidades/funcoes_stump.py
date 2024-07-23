import numpy as np


def stumpC(z):
    """
    Calcula a função Stumpff C(z).

    Parâmetros:
        z (float): Valor de entrada para a função Stumpff C.

    Retorna:
        float: Resultado da função Stumpff C(z).
    """
    # Evaluate the Stumpff function C(z)
    if z > 0:
        c = (1 - np.cos(np.sqrt(z))) / z
    elif z < 0:
        c = (np.cosh(np.sqrt(-z)) - 1) / (-z)
    else:
        c = 1 / 2

    return c


def stumpS(z):
    """
    Calcula a função Stumpff S(z).

    Parâmetros:
        z (float): Valor de entrada para a função Stumpff S.

    Retorna:
        float: Resultado da função Stumpff S(z).
    """
    if z > 0:
        return (np.sqrt(z) - np.sin(np.sqrt(z))) / (np.sqrt(z)**3)
    elif z < 0:
        return (np.sinh(np.sqrt(-z)) - np.sqrt(-z)) / (np.sqrt(-z)**3)
    else:
        return 1 / 6
