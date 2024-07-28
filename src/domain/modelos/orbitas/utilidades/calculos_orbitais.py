import numpy as np


def calcula_velocidade_orbital(parametro_gravitacional, distancia_radial, semi_eixo_maior=None):
    """
    Calcula a velocidade orbital.

    Retorna:
    float: Velocidade orbital.
    """
    if semi_eixo_maior is None:
        return np.sqrt(parametro_gravitacional / distancia_radial)
    return np.sqrt(parametro_gravitacional * ((2 / distancia_radial) - (1 / semi_eixo_maior)))


def semi_eixo_maior(periapsis, apoapsis):
    """
    Calcula o semi-eixo maior de uma órbita eliptica.

    Retorna:
    float: Semi-eixo maior.
    """
    return (periapsis + apoapsis) / 2


def excentricidade_orbital(periapsis, apoapsis):
    """
    Calcula a excentricidade orbital.

    Retorna:
    float: Excentricidade orbital.
    """
    eixo = semi_eixo_maior(periapsis, apoapsis)
    return (apoapsis / eixo) - 1


def calcula_inclinacao_heliosincrona(a, e, J2, R_e, mu, OMpm):
    # Calcular o movimento médio
    n = np.sqrt(mu / a ** 3)

    # Calcular o parâmetro orbital específico p
    p = a * (1 - e ** 2)

    # Calcular o cosseno da inclinação necessária
    cos_i = -OMpm * (2 / (3 * n * J2)) * (p / R_e) ** 2

    # Verificar se o valor calculado está dentro do domínio do arco cosseno [-1, 1]
    if cos_i < -1 or cos_i > 1:
        raise ValueError("Valor de cosseno fora do intervalo [-1, 1]. Verifique os parâmetros de entrada.")

    # Calcular a inclinação a partir do cosseno
    inclination = np.arccos(cos_i)

    # Converter a inclinação de radianos para graus
    return np.rad2deg(inclination)


def calcula_distancia_radial(posicao_celeste):
    """
    Calcula a distância radial de uma posição celeste.

    Retorna:
    float: Distância radial.
    """
    return np.linalg.norm(posicao_celeste)


def calcula_quantidade_movimento_angular_especifica(posicao_celeste, velocidade_celeste):
    """
    Calcula o momento angular.

    Retorna:
    array: Momento angular.
    """
    return np.cross(posicao_celeste, velocidade_celeste)


def calcula_excentricidade(velocidade_celeste, momento_angular, parametro_gravitacional, posicao_celeste,
                           distancia_radial):
    """
    Calcula a excentricidade da órbita.
    Retorna:
    array: Excentricidade da órbita.
    """
    return np.cross(velocidade_celeste, momento_angular) / parametro_gravitacional - posicao_celeste / distancia_radial


def calcula_parametro_e_semi_eixo_maior(momento_angular, parametro_gravitacional, excentricidade):
    """
    Calcula os parâmetros orbitais.

    Parâmetros:
    momento_angular (float): Momento angular.
    parametro_gravitacional (float): Parâmetro gravitacional.
    excentricidade (float): Excentricidade orbital.

    Retorna:
    tuple: Parâmetro e semi-eixo maior da órbita.
    """
    parametro = momento_angular ** 2 / parametro_gravitacional
    eixo = parametro / (1 - excentricidade ** 2)
    return parametro, eixo


def calcula_anomalia_verdadeira(posicao_celeste, posicao_perifocal, distancia_radial, parametro, excentricidade):
    """
    Calcula a anomalia verdadeira.

    Retorna:
    float: Anomalia verdadeira.
    """
    cos_theta = (parametro - distancia_radial) / (excentricidade * distancia_radial)
    sin_theta = np.dot(posicao_celeste, posicao_perifocal) / (distancia_radial * parametro)
    return np.arctan2(sin_theta, cos_theta)


def calcula_tempo_periastro(epoca, excentricidade, anomalia_verdadeira, parametro_gravitacional, semi_eixo_maior):
    """
    Calcula o tempo de periastro.

    Retorna:
    float: Tempo de periastro.
    """
    if excentricidade < 1:
        anomalia_excentrica = 2 * np.arctan(
            np.sqrt((1 - excentricidade) / (1 + excentricidade)) * np.tan(anomalia_verdadeira / 2))
        return epoca - (anomalia_excentrica - excentricidade * np.sin(anomalia_excentrica)) / np.sqrt(
            parametro_gravitacional / semi_eixo_maior ** 3)
    elif excentricidade == 1:
        parametro = semi_eixo_maior * (1 - excentricidade ** 2)  # Ajuste para órbita parabólica
        return epoca - (2 / 3) * ((parametro ** (3 / 2)) / np.sqrt(parametro_gravitacional)) * np.tan(
            anomalia_verdadeira / 2)
    else:
        anomalia_hiperbolica = 2 * np.arctanh(
            np.sqrt((excentricidade - 1) / (1 + excentricidade)) * np.tan(anomalia_verdadeira / 2))
        return epoca - (excentricidade * np.sinh(anomalia_hiperbolica) - anomalia_hiperbolica) / np.sqrt(
            -parametro_gravitacional / semi_eixo_maior ** 3)


def calcula_tempo_de_periastro_anomalia_verdadeira(momento_angular, excentricidade, distancia_radial, posicao_celeste,
                                                   parametro, exc, epoca, parametro_gravitacional, semi_eixo_maior):
    """
    Calcula os elementos orbitais.

    Retorna:
    tuple: Anomalia verdadeira e tempo de periastro.
    """
    posicao_perifocal = parametro * np.cross(momento_angular, excentricidade) / (np.linalg.norm(momento_angular) * exc)
    anomalia_verdadeira = calcula_anomalia_verdadeira(posicao_celeste, posicao_perifocal, distancia_radial, parametro,
                                                      exc)
    tempo_periastro = calcula_tempo_periastro(epoca, exc, anomalia_verdadeira, parametro_gravitacional, semi_eixo_maior)
    return anomalia_verdadeira, tempo_periastro


def calcula_quantidade_movimento_angular(paramentro, mu):
    return np.sqrt(paramentro * mu)


def calcula_inclinacao_raan_argumento_de_periastro(momento_angular, excentricidade):
    """
    Calcula a inclinação, RAAN e argumento de periastro.

    Retorna:
    tuple: RAAN, inclinação e argumento de periastro.
    """
    momento_angular_normalizado = momento_angular / np.linalg.norm(momento_angular)
    vetor_referencia = np.array([0, 0, 1])
    linha_dos_nodos = np.cross(vetor_referencia, momento_angular_normalizado) / np.linalg.norm(
        np.cross(vetor_referencia, momento_angular_normalizado))
    raan = np.arctan2(linha_dos_nodos[1], linha_dos_nodos[0])
    inclinacao = np.arccos(np.dot(momento_angular_normalizado, vetor_referencia))
    excentricidade_normalizada = excentricidade / np.linalg.norm(excentricidade)
    cos_argumento_de_periastro = np.dot(excentricidade_normalizada, linha_dos_nodos)
    sin_argumento_de_periastro = np.dot(momento_angular_normalizado,
                                        np.cross(linha_dos_nodos, excentricidade_normalizada))
    argumento_de_periastro = np.arctan2(sin_argumento_de_periastro, cos_argumento_de_periastro)
    return raan, inclinacao, argumento_de_periastro


def determina_parametros_orbitais(tempo_observacao, parametro_gravitacional, posicao_celeste, velocidade_celeste):
    """
    Determina os parâmetros orbitais.

    Algoritmos baseados em: Tewari, Ashish. (2011). Atmospheric and space flight dynamics

    Retorna:
    objeto: Órbita.
    """
    posicao_celeste = np.asarray(posicao_celeste, dtype=np.float64).flatten()
    velocidade_celeste = np.asarray(velocidade_celeste, dtype=np.float64).flatten()
    distancia_radial = calcula_distancia_radial(posicao_celeste)
    momento_angular = calcula_quantidade_movimento_angular_especifica(posicao_celeste, velocidade_celeste)
    excentricidade = calcula_excentricidade(velocidade_celeste, momento_angular, parametro_gravitacional,
                                            posicao_celeste, distancia_radial)
    exc = np.linalg.norm(excentricidade)
    h = np.linalg.norm(momento_angular)
    parametro, eixo = calcula_parametro_e_semi_eixo_maior(h, parametro_gravitacional, exc)
    anomalia_verdadeira, tempo_periastro = calcula_tempo_de_periastro_anomalia_verdadeira(momento_angular,
                                                                                          excentricidade,
                                                                                          distancia_radial,
                                                                                          posicao_celeste, parametro,
                                                                                          exc, tempo_observacao,
                                                                                          parametro_gravitacional, eixo)
    raan, inclinacao, argumento_de_periastro = calcula_inclinacao_raan_argumento_de_periastro(momento_angular,
                                                                                              excentricidade)

    return eixo, exc, inclinacao, raan, argumento_de_periastro, anomalia_verdadeira, tempo_periastro


def calcular_periodo_orbital(a, mu):
    """
    Calcula o período orbital.

    Parameters:
    semi_eixo_maior (float): O semi-eixo maior da órbita (km).
    mu (float): O parâmetro gravitacional (km^3/s^2).

    Returns:
    float: O período orbital (s).
    """
    return 2 * np.pi / np.sqrt(mu) * a**1.5
