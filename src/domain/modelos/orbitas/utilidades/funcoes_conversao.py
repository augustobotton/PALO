import numpy as np


def converte_longitude_fixa_planeta_para_longitude_celeste(
        tempo_longitude_celeste_desejada: float,
        longitude_fixa_planeta: float,
        velocidade_rotacao_planeta: float,
        tempo_meridiano_ref_long_celeste_nula: float
) -> float:
    """
    Calcula a longitude celeste a partir da longitude fixa ao planeta.
    :return: Longitude celeste no tempo desejado (rad).
    """
    longitude_relativa_ref_fixo_planeta = longitude_fixa_planeta + velocidade_rotacao_planeta * (
            tempo_longitude_celeste_desejada - tempo_meridiano_ref_long_celeste_nula)
    return longitude_relativa_ref_fixo_planeta


def rvel_polar_para_rvel_retangular(
        modulo_vetor_velocidade: float,
        azimute_da_velocidade: float,
        elevacao_da_velocidade: float,
        distancia_radial: float,
        latitude: float,
        longitude: float
) -> tuple[np.ndarray, np.ndarray]:
    """
    Converte a velocidade do sistema LVLH (coordenadas polares) para o sistema ECI ou ECEF retangular.

    :param modulo_vetor_velocidade: Módulo do vetor velocidade (m/s).
    :param azimute_da_velocidade: Azimute da velocidade (rad).
    :param elevacao_da_velocidade: Elevação da velocidade (rad).
    :param distancia_radial: Distância radial (m).
    :param latitude: Latitude (rad).
    :param longitude: Longitude no referencial desejado (ECI ou ECEF) (rad).
    :return: Vetores posição e velocidade em coordenadas retangulares no sistema ECI ou ECEF.
    """
    matriz_converte_para_LVLH = np.array([
        [np.cos(latitude) * np.cos(longitude), np.cos(latitude) * np.sin(longitude), np.sin(latitude)],
        [-np.sin(longitude), np.cos(longitude), 0],
        [-np.sin(latitude) * np.cos(longitude), -np.sin(latitude) * np.sin(longitude), np.cos(latitude)]
    ], dtype=('object'))

    vetor_velocidade_cood_cartezianas_lvlh = modulo_vetor_velocidade * np.array([
        [np.sin(elevacao_da_velocidade)],
        [np.cos(elevacao_da_velocidade) * np.sin(azimute_da_velocidade)],
        [np.cos(elevacao_da_velocidade) * np.cos(azimute_da_velocidade)]
    ], dtype=('object'))

    vetor_velocidade_coord_retangular = matriz_converte_para_LVLH.T @ vetor_velocidade_cood_cartezianas_lvlh

    vetor_posicao_LVLH = np.array([[distancia_radial], [0], [0]])
    vetor_posicao_coord_retangular = matriz_converte_para_LVLH.T @ vetor_posicao_LVLH

    return vetor_posicao_coord_retangular, vetor_velocidade_coord_retangular


def componentes_vel_relativa_para_inercial(velocidade_relativa, inclinacao_velocidade_relativa, azimute_da_velocidade_relativa,
                                 velocidade_rotacao_ref_girante, distancia_radial_origem_inercial, latitude):
    # Cálculos do ângulo de azimute da velocidade inercial
    angulo_azimute_velocidade_inercial = np.arctan2(
        velocidade_relativa * np.cos(inclinacao_velocidade_relativa) * np.sin(azimute_da_velocidade_relativa) +
        velocidade_rotacao_ref_girante * distancia_radial_origem_inercial * np.cos(latitude),
        velocidade_relativa * np.cos(inclinacao_velocidade_relativa) * np.cos(azimute_da_velocidade_relativa)
    )
    if angulo_azimute_velocidade_inercial < 0:
        angulo_azimute_velocidade_inercial += 2 * np.pi

    # Cálculo da magnitude da velocidade referencial inercial
    magnitude_velocidade_referencial_inercial = np.sqrt(
        velocidade_relativa ** 2 +
        2 * velocidade_relativa * np.cos(inclinacao_velocidade_relativa) * np.sin(azimute_da_velocidade_relativa) *
        distancia_radial_origem_inercial * velocidade_rotacao_ref_girante * np.cos(latitude) +
        distancia_radial_origem_inercial ** 2 * velocidade_rotacao_ref_girante ** 2 * np.cos(latitude) ** 2
    )

    # Cálculo do ângulo de trajetória
    angulo_de_trajetoria = np.arctan2(
        np.sin(inclinacao_velocidade_relativa) * np.cos(angulo_azimute_velocidade_inercial),
        np.cos(inclinacao_velocidade_relativa) * np.cos(azimute_da_velocidade_relativa)
    )
    if abs(angulo_de_trajetoria) > np.pi / 2:
        if angulo_de_trajetoria < np.pi / 2:
            angulo_de_trajetoria += np.pi
        if angulo_de_trajetoria > np.pi / 2:
            angulo_de_trajetoria -= np.pi

    return magnitude_velocidade_referencial_inercial, angulo_de_trajetoria, angulo_azimute_velocidade_inercial


def matriz_rotacao_orbital_inercial(raan: float, arg_perigeu: float, inclinacao: float) -> np.ndarray:
    """
    Cria a matriz de rotação do plano orbital (perifocal) para o sistema de coordenadas inerciais (celeste).

    :return: Matriz de rotação 3x3.
    """
    cos_raan = np.cos(raan)
    sin_raan = np.sin(raan)
    cos_arg_perigeu = np.cos(arg_perigeu)
    sin_arg_perigeu = np.sin(arg_perigeu)
    cos_inclinacao = np.cos(inclinacao)
    sin_inclinacao = np.sin(inclinacao)

    matriz_rotacao = np.array([
        [
            cos_raan * cos_arg_perigeu - sin_raan * sin_arg_perigeu * cos_inclinacao,
            -cos_raan * sin_arg_perigeu - sin_raan * cos_arg_perigeu * cos_inclinacao,
            sin_raan * sin_inclinacao
        ],
        [
            sin_raan * cos_arg_perigeu + cos_raan * sin_arg_perigeu * cos_inclinacao,
            -sin_raan * sin_arg_perigeu + cos_raan * cos_arg_perigeu * cos_inclinacao,
            -cos_raan * sin_inclinacao
        ],
        [
            sin_arg_perigeu * sin_inclinacao,
            cos_arg_perigeu * sin_inclinacao,
            cos_inclinacao
        ]
    ])

    return matriz_rotacao


# Matrizes de rotação elementares
def rotacao_x(angulo: float) -> np.ndarray:
    """
    Gera a matriz de rotação ao redor do eixo X por um ângulo.

    :param angulo: Ângulo de rotação em radianos.
    :return: Matriz de rotação 3x3.
    """
    matriz_rotacao = np.array([
        [1, 0, 0],
        [0, np.cos(angulo), np.sin(angulo)],
        [0, -np.sin(angulo), np.cos(angulo)]
    ])
    return matriz_rotacao


def rotacao_y(angulo: float) -> np.ndarray:
    """
    Gera a matriz de rotação ao redor do eixo Y por um ângulo.

    :param angulo: Ângulo de rotação em radianos.
    :return: Matriz de rotação 3x3.
    """
    matriz_rotacao = np.array([
        [np.cos(angulo), 0, -np.sin(angulo)],
        [0, 1, 0],
        [np.sin(angulo), 0, np.cos(angulo)]
    ])
    return matriz_rotacao


def rotacao_z(angulo: float) -> np.ndarray:
    """
    Gera a matriz de rotação ao redor do eixo Z por um ângulo.

    :param angulo: Ângulo de rotação em radianos.
    :return: Matriz de rotação 3x3.
    """
    matriz_rotacao = np.array([
        [np.cos(angulo), np.sin(angulo), 0],
        [-np.sin(angulo), np.cos(angulo), 0],
        [0, 0, 1]
    ])
    return matriz_rotacao


def matriz_inercial_para_perifocal(
        longitude_no_ascendente: float,
        inclinacao_orbita: float,
        argumento_perigeu: float
) -> np.ndarray:
    """
    Calcula a matriz de rotação do sistema de referência inercial para o sistema perifocal.

    :param longitude_no_ascendente: Longitude do nó ascendente em radianos.
    :param inclinacao_orbita: Inclinação da órbita em radianos.
    :param argumento_perigeu: Argumento do perigeu em radianos.
    :return: Matriz de rotação 3x3 do sistema inercial para o perifocal.
    """
    matriz_rotacao = rotacao_z(argumento_perigeu) @ rotacao_x(inclinacao_orbita) @ rotacao_z(longitude_no_ascendente)
    return matriz_rotacao


def perifocal_para_inercial(vector, OMEGA, i, omega):
    """
    Converte um vetor do referencial perifocal para o referencial inercial.

    Parâmetros:
    vector (np.array): Vetor no referencial perifocal.
    OMEGA (float): Longitude do nodo ascendente em radianos.
    i (float): Inclinação da órbita em radianos.
    omega (float): Argumento do periastro em radianos.

    Retorna:
    np.array: Vetor no referencial inercial.
    """
    # Matriz de rotação composta
    matriz_rotacao = np.array([
        [np.cos(OMEGA) * np.cos(omega) - np.sin(OMEGA) * np.sin(omega) * np.cos(i),
         -np.cos(OMEGA) * np.sin(omega) - np.sin(OMEGA) * np.cos(omega) * np.cos(i),
         np.sin(OMEGA) * np.sin(i)],
        [np.sin(OMEGA) * np.cos(omega) + np.cos(OMEGA) * np.sin(omega) * np.cos(i),
         -np.sin(OMEGA) * np.sin(omega) + np.cos(OMEGA) * np.cos(omega) * np.cos(i),
         -np.cos(OMEGA) * np.sin(i)],
        [np.sin(omega) * np.sin(i),
         np.cos(omega) * np.sin(i),
         np.cos(i)]
    ])

    # Aplica a matriz de rotação ao vetor
    return np.dot(matriz_rotacao, vector)
