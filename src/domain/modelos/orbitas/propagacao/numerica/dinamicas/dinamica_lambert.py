import numpy as np
def dinamica_lambert(t, f, *args):
    """
    Calcula a dinâmica de um problema de Lambert.

    :param t: Tempo (s).
    :param f: Vetor de estado, 6x1.
    :param args: Parâmetro gravitacional (m^3/s^2).
    :return: Lista com as componentes da velocidade e da aceleração.
    """
    mu = args[0]
    x, y, z, vx, vy, vz = f
    r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
    ax = -mu * x / r ** 3
    ay = -mu * y / r ** 3
    az = -mu * z / r ** 3
    return [vx, vy, vz, ax, ay, az]
