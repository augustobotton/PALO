import numpy as np


def converte_longitudeFixaPlaneta_para_longitude_celeste(tempo_longitude_celeste_desejada,
                                                         longitude_fixaPlaneta,
                                                         velocidade_rotacao_planeta,
                                                         tempo_meridianoRef_longCeleste_nula):
    # Função para calcular a longitude celeste a partir da longitude fixa ao planeta
    # Entradas
    # t (s) - Tempo no qual se deseja saber a longitude celeste
    # long (rad) - Longitude relativa ao referencial fixo ao planeta
    # we (rad/s) - Velocidade de rotação do planeta
    # tg (s) - Tempo no qual o meridiano de referência tem longitude celeste nula
    # Saída
    # long_c (rad) - Longitude celeste no tempo t

    longitude_relativa_ref_fixo_planeta = longitude_fixaPlaneta + velocidade_rotacao_planeta * (
            tempo_longitude_celeste_desejada - tempo_meridianoRef_longCeleste_nula)
    return longitude_relativa_ref_fixo_planeta


def RvelPolar2RvelRet(modulo_vetor_velocidade, azimute_da_velocidade, elevacao_da_velocidade,
                      distancia_radial, latitude, longitude):
    """
        Função para converter velocidade do sistema LVLH (coordenadas polares) para o sistema ECI ou ECEF
        retangular
        """
    # Entradas:
    # v (m/s): Módulo do vetor velocidade
    # A (rad): Azimute da velocidade
    # phi (rad): Elevação da velocidade
    # r (m): distância radial
    # lat (rad): Latitude
    # lon (rad): Longitude no referencial desejado (ECI ou ECEF)
    # Saídas:
    # R=[R_X, R_Y, R_Z]^T (m): Vetor posição em coordenadas retangulares no sistema ECI ou ECEF (dependendo
    # da entrada de dados de longitude)
    # V=[V_X, V_Y, V_Z]^T (m/s): Vetor velocidade em coordenadas retangulares no sistema ECI ou ECEF (
    # dependendo da entrada de dados de longitude). Pode ser a velocidade relativa ou a inercial,
    # dependendo dos dados de velocidade fornecidos.

    ## Cálculos
    # Matriz de conversão do sistema ECI ou ECEF para o LVLH

    matriz_converte_para_LVLH = np.array(
        [[np.cos(latitude) * np.cos(longitude), np.cos(latitude) * np.sin(longitude), np.sin(latitude)],
         [-np.sin(longitude), np.cos(longitude), 0],
         [-np.sin(latitude) * np.cos(longitude), -np.sin(latitude) * np.sin(longitude),
          np.cos(latitude)]], dtype=('object'));
    # Vetor velocidade em coordenadas cartezianas no sistema LVLH
    vetor_velocidade_cood_cartezianas_lvlh = modulo_vetor_velocidade * np.array(
        [[np.sin(elevacao_da_velocidade)],
         [np.cos(elevacao_da_velocidade) * np.sin(azimute_da_velocidade)],
         [np.cos(elevacao_da_velocidade) * np.cos(azimute_da_velocidade)]], dtype=('object'));
    # Transformação da velocidade para o sistema ECI ou ECEF em coordenadas retangulares

    vetor_velocidade_coord_retangular = matriz_converte_para_LVLH.T @ vetor_velocidade_cood_cartezianas_lvlh[
                                                                      :][:]

    # Vetor posição no sistema LHVLH
    vetor_posicao_LVLH = np.array([[distancia_radial],
                                   [0],
                                   [0]])
    # Transformação da posição para o sistema ECI ou ECEF em coordenadas
    # retangulares
    vetor_posicao_coord_retangular = matriz_converte_para_LVLH.T @ vetor_posicao_LVLH
    return vetor_posicao_coord_retangular, vetor_velocidade_coord_retangular


def Vrel2Vine(velocidade_relativa, inclinacao_velocidade_relativa, azimute_da_velocidade_relativa,
              velocidade_rotacao_ref_girante, distancia_radial_origem_inercial, latitude):
    # Função para converter a velocidade, elevação e azimute da velocidade relativa para a velocidade
    # inercial
    # Entradas
    # vr (m/s): velocidade relativa
    # phir (rad): inclinação da velocidade relativa
    # Ar (rad): azimute da velocidade relativa
    # we (rad/s): velocidade de rotação do referencial girante
    # r (m): distância radial até a origem do referencial inercial
    # dt (rad): latitude
    # Saídas
    # v (m/s): magnitude da velocidade com respeito ao referencial inercial
    # phi (rad): ângulo de elevação da velocidade inercial (angulo de trajetoria)
    # A (rad): ângulo de azimute da velocidade inercial
    ## Cálculos

    angulo_azimute_velocidade_inercial = np.arctan2(velocidade_relativa * np.cos(
        inclinacao_velocidade_relativa) * np.sin(
        azimute_da_velocidade_relativa) + velocidade_rotacao_ref_girante * distancia_radial_origem_inercial
                                                    * np.cos(
        latitude), velocidade_relativa * np.cos(inclinacao_velocidade_relativa) * np.cos(
        azimute_da_velocidade_relativa))
    if angulo_azimute_velocidade_inercial < 0:
        angulo_azimute_velocidade_inercial += 2 * np.pi
    magnitude_velocidade_referencial_inercial = np.sqrt(
        velocidade_relativa ** 2 + 2 * velocidade_relativa * np.cos(inclinacao_velocidade_relativa) * np.sin(
            azimute_da_velocidade_relativa) * distancia_radial_origem_inercial *
        velocidade_rotacao_ref_girante * np.cos(
            latitude) + distancia_radial_origem_inercial ** 2 * velocidade_rotacao_ref_girante ** 2 * np.cos(
            latitude) ** 2)
    angulo_de_trajetoria = np.arctan(
        (np.sin(inclinacao_velocidade_relativa) * np.cos(angulo_azimute_velocidade_inercial)) / (
                np.cos(inclinacao_velocidade_relativa) * np.cos(azimute_da_velocidade_relativa)))

    return (magnitude_velocidade_referencial_inercial, angulo_de_trajetoria,
            angulo_azimute_velocidade_inercial)


def matriz_rotacao_orbital_inercial(raan: float, arg_perigeu: float, inclinacao: float) -> np.ndarray:
    """
    Cria a matriz de rotação do plano orbital (perifocal) para o sistema de coordenadas inerciais (
    celeste).

    :param raan: Longitude do nó ascendente (em radianos).
    :param arg_perigeu: Argumento do perigeu (em radianos).
    :param inclinacao: Inclinação da órbita (em radianos).
    :return: Matriz de rotação 3x3.
    """
    cos_raan = np.cos(raan)
    sin_raan = np.sin(raan)
    cos_arg_perigeu = np.cos(arg_perigeu)
    sin_arg_perigeu = np.sin(arg_perigeu)
    cos_inclinacao = np.cos(inclinacao)
    sin_inclinacao = np.sin(inclinacao)

    matriz_rotacao = np.array([
        [cos_raan * cos_arg_perigeu - sin_raan * sin_arg_perigeu * cos_inclinacao,
         -cos_raan * sin_arg_perigeu - sin_raan * cos_arg_perigeu * cos_inclinacao,
         sin_raan * sin_inclinacao],
        [sin_raan * cos_arg_perigeu + cos_raan * sin_arg_perigeu * cos_inclinacao,
         -sin_raan * sin_arg_perigeu + cos_raan * cos_arg_perigeu * cos_inclinacao,
         -cos_raan * sin_inclinacao],
        [sin_arg_perigeu * sin_inclinacao,
         cos_arg_perigeu * sin_inclinacao,
         cos_inclinacao]
    ])

    return matriz_rotacao
