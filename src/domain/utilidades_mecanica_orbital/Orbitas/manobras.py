import numpy as np

from src.domain.utilidades_mecanica_orbital.Orbitas.ModeloOrbita import Orbita
from src.domain.utilidades_mecanica_orbital.Orbitas.orbitalUtils import Converte
from src.domain.utilidades_mecanica_orbital.Orbitas.orbitalUtils.calculos_orbitais import *

def manobra_hohmann(orbita_inicial, raio_da_orbita_desejada):

    r_i = np.linalg.norm(orbita_inicial.calcular_vetor_posicao_orbital())
    v_i = np.linalg.norm(orbita_inicial.calcular_vetor_velocidade_orbital())
    a_trans = (r_i + raio_da_orbita_desejada) / 2

    dv_a = np.sqrt(2 * orbita_inicial.mu / r_i - orbita_inicial.mu / a_trans) - v_i
    dv_b = np.sqrt(orbita_inicial.mu / raio_da_orbita_desejada) - np.sqrt(2 * orbita_inicial.mu / raio_da_orbita_desejada - orbita_inicial.mu / a_trans)

    dv_a = np.array([0, dv_a, 0])
    dv_b = np.array([0, -dv_b, 0])

    matriz_de_rotacao = Converte.matriz_rotacao_orbital_inercial(orbita_inicial.raan, orbita_inicial.arg_periastro, orbita_inicial.inclinacao)

    dv_a = matriz_de_rotacao @ dv_a
    dv_b = matriz_de_rotacao @ dv_b

    t_trans = np.pi * np.sqrt(a_trans ** 3 / orbita_inicial.mu)

    return dv_a, dv_b, t_trans

def manobra_mono_impulsiva(orbita_inicial: Orbita, orbita_final: Orbita):

     semi_eixo_orbita_final = semi_eixo_maior(orbita_final.calcula_periastro(), orbita_final.calcula_apoastro())
     r_f = np.linalg.norm(orbita_inicial.calcular_vetor_posicao_orbital())
     v_f = calcaula_velociade_orbital(orbita_final.mu,r_f,semi_eixo_orbita_final)
     h_f = calcula_quantidade_movimento_angular()


     return v_f









def bielliptic(k, r_b, r_f, rv):
    r"""Calculate the increments in the velocities and the time of flight of the maneuver.

    The bielliptic maneuver employs two Hohmann transfers, therefore two
    intermediate orbits are established. We define the different radius
    relationships as follows:

    .. math::
        \begin{align}
            a_{trans1} &= \frac{r_{i} + r_{b}}{2}\\
            a_{trans2} &= \frac{r_{b} + r_{f}}{2}\\
        \end{align}

    The increments in the velocity are:

    .. math::
        \begin{align}
            \Delta v_{a} &= \sqrt{\frac{2\mu}{r_{i}} - \frac{\mu}{a_{trans1}}} - v_{i}\\
            \Delta v_{b} &= \sqrt{\frac{2\mu}{r_{b}} - \frac{\mu}{a_{trans2}}} - \sqrt{\frac{2\mu}{r_{b}} - \frac{\mu}{a_trans{1}}}\\
            \Delta v_{c} &= \sqrt{\frac{\mu}{r_{f}}} - \sqrt{\frac{2\mu}{r_{f}} - \frac{\mu}{a_{trans2}}}\\
        \end{align}

    The time of flight for this maneuver is the addition of the time needed for both transition orbits, following the same formula as
    Hohmann:

    .. math::
        \begin{align}
            \tau_{trans1} &= \pi \sqrt{\frac{a_{trans1}^{3}}{\mu}}\\
            \tau_{trans2} &= \pi \sqrt{\frac{a_{trans2}^{3}}{\mu}}\\
        \end{align}

    Parameters
    ----------
    k : float
        Standard Gravitational parameter
    r_b : float
        Altitude of the intermediate orbit
    r_f : float
        Final orbital radius
    rv : numpy.ndarray, numpy.ndarray
        Position and velocity vectors

    """
    _, ecc, inc, raan, argp, nu = rv2coe(k, *rv)
    h_i = norm(cross(*rv))
    p_i = h_i**2 / k

    r_i, v_i = rv_pqw(k, p_i, ecc, nu)

    r_i = norm(r_i)
    v_i = norm(v_i)
    a_trans1 = (r_i + r_b) / 2
    a_trans2 = (r_b + r_f) / 2

    dv_a = np.sqrt(2 * k / r_i - k / a_trans1) - v_i
    dv_b = np.sqrt(2 * k / r_b - k / a_trans2) - np.sqrt(
        2 * k / r_b - k / a_trans1
    )
    dv_c = np.sqrt(k / r_f) - np.sqrt(2 * k / r_f - k / a_trans2)

    dv_a = np.array([0, dv_a, 0])
    dv_b = np.array([0, -dv_b, 0])
    dv_c = np.array([0, dv_c, 0])

    rot_matrix = coe_rotation_matrix(inc, raan, argp)

    dv_a = rot_matrix @ dv_a
    dv_b = rot_matrix @ dv_b
    dv_c = rot_matrix @ dv_c

    t_trans1 = np.pi * np.sqrt(a_trans1**3 / k)
    t_trans2 = np.pi * np.sqrt(a_trans2**3 / k)

    return dv_a, dv_b, dv_c, t_trans1, t_trans2