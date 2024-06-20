from src.domain.utilidades_mecanica_orbital.Orbitas.ModeloOrbita import Orbita
from src.domain.utilidades_mecanica_orbital.orbitalUtils import Converte
from src.domain.utilidades_mecanica_orbital.orbitalUtils.calculos_orbitais import *


def manobra_hohmann(orbita_inicial, raio_da_orbita_desejada):
    r_i = np.linalg.norm(orbita_inicial.calcular_vetor_posicao_orbital())
    v_i = np.linalg.norm(orbita_inicial.calcular_vetor_velocidade_orbital())
    a_trans = (r_i + raio_da_orbita_desejada) / 2

    dv_a = np.sqrt(2 * orbita_inicial.mu / r_i - orbita_inicial.mu / a_trans) - v_i
    dv_b = np.sqrt(orbita_inicial.mu / raio_da_orbita_desejada) - np.sqrt(
        2 * orbita_inicial.mu / raio_da_orbita_desejada - orbita_inicial.mu / a_trans)

    dv_a = np.array([0, dv_a, 0])
    dv_b = np.array([0, -dv_b, 0])

    matriz_de_rotacao = Converte.matriz_rotacao_orbital_inercial(orbita_inicial.raan, orbita_inicial.arg_periastro,
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

    matriz_de_rotacao = Converte.matriz_rotacao_orbital_inercial(orbita_inicial.raan, orbita_inicial.arg_periastro,
                                                                 orbita_inicial.inclinacao)

    dv_a = matriz_de_rotacao @ dv_a
    dv_b = matriz_de_rotacao @ dv_b
    dv_c = matriz_de_rotacao @ dv_c

    t_trans1 = np.pi * np.sqrt(a_trans1 ** 3 / orbita_inicial.mu)
    t_trans2 = np.pi * np.sqrt(a_trans2 ** 3 / orbita_inicial.mu)

    return dv_a, dv_b, dv_c, t_trans1, t_trans2


def manobra_coplanar(orbita_inicial: Orbita, orbita_final: Orbita):
    semi_eixo_orbita_final = semi_eixo_maior(orbita_final.calcula_periastro(), orbita_final.calcula_apoastro())
    r_f = np.linalg.norm(orbita_inicial.calcular_vetor_posicao_orbital())
    v_i = calcaula_velocidade_orbital(orbita_inicial.mu, r_f)
    v_f = calcaula_velocidade_orbital(orbita_final.mu, r_f, semi_eixo_orbita_final)

    # Considerar possibilidade de orbita inicial eliptica
    if orbita_inicial.excentricidade > 0:

        anomalia_verdadeira = calcula_anomalia_verdadeira(orbita_inicial, r_f) #TODO conferir equação para anomalia verdadeira
        angulo_alpha = np.arccos(
            np.cos(anomalia_verdadeira) / (1 + orbita_inicial.excentricidade * np.cos(anomalia_verdadeira)))
    else:
        angulo_alpha = np.arccos((r_f * v_i) / (r_f * v_f))

    delta_v = delta_velocidade(v_i, v_f, angulo_alpha)
    betha = np.arcsin(v_f / delta_v * np.sin(angulo_alpha))

    return delta_v, betha


def manobra_mudanca_de_plano(delta_inclinacao, vi):
    delta_v = 2 * vi * np.sin(delta_inclinacao / 2)
    betha = np.pi / 2 + delta_inclinacao / 2

    return delta_v, betha


def manobra_mudanca_de_plano_generica():
    return 0