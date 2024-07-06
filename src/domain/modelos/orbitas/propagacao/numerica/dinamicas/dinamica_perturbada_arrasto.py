import numpy as np

from src.domain.modelos.orbitas.utilidades.densidade_altitude import densidade_altitude


def dinamica_perturbada_arrasto(t, f, *args):
    """
    Calcula as taxas de variação da posição e da velocidade.

    Argumentos:
        t  - Tempo (s)
        f  - Vetor de estado [x, y, z, vx, vy, vz] (km e km/s)
        RE - Raio da Terra (km)
        wE - Vetor de velocidade angular da Terra (rad/s)
        CD - Coeficiente de arrasto
        A  - Área transversal do objeto (m^2)
        m  - Massa do objeto (kg)
        mu - Parâmetro gravitacional (km^3/s^2)

    Retorna:
        dfdt - Vetor de derivadas [vx, vy, vz, ax, ay, az] (km/s e km/s^2)
    """

    planeta, m, A, CD = args
    wE = np.array([0, 0, 7.2921159e-5])  # Velocidade angular da Terra (rad/s)

    # Vetor posição (km)
    R = f[:3]

    r = np.linalg.norm(R)  # Distância ao centro da Terra (km)
    alt = r - planeta.raio_equatorial  # Altitude (km)

    # Vetor velocidade (km/s)
    V = f[3:]

    Vrel = V - np.cross(wE, R)  # Velocidade relativa à atmosfera (km/s)
    vrel = np.linalg.norm(Vrel)  # Velocidade relativa à atmosfera (km/s)
    uv = Vrel / vrel  # Vetor unitário da velocidade relativa
    densidade_do_ar = densidade_altitude(alt)

    # Aceleração devido ao arrasto (m/s^2)
    ap = -(CD * A / m) * densidade_do_ar * (vrel * 1000) ** 2 / 2 * uv  # Convertendo vrel de km/s para m/s
    a0 = -planeta.mut * R / r ** 3  # Aceleração gravitacional (km/s^2)

    a = a0 + (ap / 1000)  # Aceleração total (km/s^2), convertendo ap de m/s^2 para km/s^2

    return np.concatenate([V, a])
