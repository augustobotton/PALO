import numpy as np

from src.domain.modelos.orbitas.utilidades.calculos_orbitais import calcula_velocidade_orbital
from src.domain.modelos.orbitas.utilidades.funcoes_conversao import matriz_rotacao_orbital_inercial


def manobra_hohmann(orbita_inicial, raio_da_orbita_desejada):
    r_i = np.linalg.norm(orbita_inicial.calcular_vetor_posicao_orbital())
    v_i = np.linalg.norm(orbita_inicial.calcular_vetor_velocidade_orbital())
    a_trans = (r_i + raio_da_orbita_desejada) / 2

    dv_a = np.sqrt(2 * orbita_inicial.mu / r_i - orbita_inicial.mu / a_trans) - v_i
    dv_b = np.sqrt(orbita_inicial.mu / raio_da_orbita_desejada) - np.sqrt(
        2 * orbita_inicial.mu / raio_da_orbita_desejada - orbita_inicial.mu / a_trans)

    dv_a = np.array([0, dv_a, 0])
    dv_b = np.array([0, -dv_b, 0])

    matriz_de_rotacao = matriz_rotacao_orbital_inercial(orbita_inicial.raan, orbita_inicial.arg_periastro,
                                                        orbita_inicial.inclinacao)

    dv_a = matriz_de_rotacao @ dv_a
    dv_b = matriz_de_rotacao @ dv_b

    t_trans = np.pi * np.sqrt(a_trans ** 3 / orbita_inicial.mu)

    return dv_a, dv_b, t_trans


def manobra_bieliptica(raio_intermediario, raio_da_orbita_desejada, orbita_inicial):
    r_i = np.linalg.norm(orbita_inicial.calcular_vetor_posicao_orbital())
    v_i = np.linalg.norm(orbita_inicial.calcular_vetor_velocidade_orbital())

    a_trans1 = (r_i + raio_intermediario) / 2
    a_trans2 = (raio_intermediario + raio_da_orbita_desejada) / 2

    dv_a = np.sqrt(2 * orbita_inicial.mu / r_i - orbita_inicial.mu / a_trans1) - v_i
    dv_b = np.sqrt(2 * orbita_inicial.mu / raio_intermediario - orbita_inicial.mu / a_trans2) - np.sqrt(
        2 * orbita_inicial.mu / raio_intermediario - orbita_inicial.mu / a_trans1
    )
    dv_c = np.sqrt(orbita_inicial.mu / raio_da_orbita_desejada) - np.sqrt(
        2 * orbita_inicial.mu / raio_da_orbita_desejada - orbita_inicial.mu / a_trans2)

    dv_a = np.array([0, dv_a, 0])
    dv_b = np.array([0, -dv_b, 0])
    dv_c = np.array([0, dv_c, 0])

    matriz_de_rotacao = matriz_rotacao_orbital_inercial(orbita_inicial.raan, orbita_inicial.arg_periastro,
                                                        orbita_inicial.inclinacao)

    dv_a = matriz_de_rotacao @ dv_a
    dv_b = matriz_de_rotacao @ dv_b
    dv_c = matriz_de_rotacao @ dv_c

    t_trans1 = np.pi * np.sqrt(a_trans1 ** 3 / orbita_inicial.mu)
    t_trans2 = np.pi * np.sqrt(a_trans2 ** 3 / orbita_inicial.mu)

    return dv_a, dv_b, dv_c, t_trans1, t_trans2


def manobra_mudanca_de_plano(delta_inclinacao, vi):
    delta_v = 2 * vi * np.sin(delta_inclinacao / 2)
    betha = np.pi / 2 + delta_inclinacao / 2

    return delta_v, betha

def delta_velocidade(velocidade_inicial, velocidade_final, angulo=None):
    """
    Calcula a variação de velocidade (delta-v).

    Retorna:
    float: Variação de velocidade (delta-v).
    """
    if angulo is None:
        return abs(velocidade_final - velocidade_inicial)
    cosphi = np.cos(np.deg2rad(angulo))
    return np.sqrt(velocidade_inicial ** 2 + velocidade_final ** 2 - 2 * velocidade_inicial * velocidade_final * cosphi)

def aplicar_delta_v(v, delta_v):
    """
    Aplica o delta_v ao vetor de velocidade atual.

    Parâmetros:
    v (ndarray): Vetor de velocidade atual.
    delta_v (float): Valor de delta_v a ser aplicado (km/s).

    Retorna:
    ndarray: Novo vetor de velocidade após aplicar delta_v.
    """
    # Considera delta_v na direção da velocidade atual
    v_unitario = v / np.linalg.norm(v)
    return v + delta_v * v_unitario




def delta_v_perigee_raise(r_p, r_a, r_p_prime, mu):
    """
    Calcula o delta_v necessário para uma manobra de aumento de perigeu.

    Parâmetros:
    r_p (float): Raio do perigeu inicial (km).
    r_a (float): Raio do apogeu (km).
    r_p_prime (float): Novo raio do perigeu desejado (km).
    mu (float): Parâmetro gravitacional do corpo central (km^3/s^2).

    Retorna:
    float: Delta_v necessário para a manobra (km/s).
    """
    # Semi-eixo maior da órbita inicial
    a = (r_p + r_a) / 2

    # Semi-eixo maior da nova órbita
    a_prime = (r_p_prime + r_a) / 2

    # Velocidade orbital no apogeu da órbita inicial
    v_a = calcula_velocidade_orbital(mu,r_a, a)

    # Velocidade orbital no apogeu da nova órbita
    v_a_prime = calcula_velocidade_orbital(mu, r_a, a_prime)

    # Delta_v necessário
    delta_v = v_a_prime - v_a

    return delta_v