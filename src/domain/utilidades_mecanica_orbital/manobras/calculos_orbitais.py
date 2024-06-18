import numpy as np

from src.domain.utilidades_mecanica_orbital.Orbitas.ConstrutorOrbita import ConstrutorOrbita

# Constantes globais
raio_da_terra = 6378.137 * 1000  # [km]
parametro_gravitacional = 3.986e14  # [m^2/s^2]

def velocidade_orbital(raio, semi_eixo_maior=None):
    """
    Calcula a velocidade orbital.

    Parâmetros:
    raio (float): Distância radial.
    semi_eixo_maior (float, opcional): Semi-eixo maior da órbita.

    Retorna:
    float: Velocidade orbital.
    """
    if semi_eixo_maior is None:
        return np.sqrt(parametro_gravitacional / raio)
    return np.sqrt(parametro_gravitacional * ((2 / raio) - (1 / semi_eixo_maior)))

def semi_eixo_maior(periapsis, apoapsis):
    """
    Calcula o semi-eixo maior de uma órbita.

    Parâmetros:
    periapsis (float): Distância do periapsis.
    apoapsis (float): Distância do apoapsis.

    Retorna:
    float: Semi-eixo maior.
    """
    return (periapsis + apoapsis) / 2

def excentricidade_orbital(periapsis, apoapsis):
    """
    Calcula a excentricidade orbital.

    Parâmetros:
    periapsis (float): Distância do periapsis.
    apoapsis (float): Distância do apoapsis.

    Retorna:
    float: Excentricidade orbital.
    """
    eixo = semi_eixo_maior(periapsis, apoapsis)
    return (apoapsis / eixo) - 1

def delta_velocidade(velocidade_inicial, velocidade_final, angulo=None):
    """
    Calcula a variação de velocidade (delta-v).

    Parâmetros:
    velocidade_inicial (float): Velocidade inicial.
    velocidade_final (float): Velocidade final.
    angulo (float, opcional): Ângulo entre as velocidades.

    Retorna:
    float: Variação de velocidade (delta-v).
    """
    if angulo is None:
        return abs(velocidade_final - velocidade_inicial)
    cosphi = np.cos(np.deg2rad(angulo))
    return np.sqrt(velocidade_inicial**2 + velocidade_final**2 - 2 * velocidade_inicial * velocidade_final * cosphi)

def calcula_distancia_radial(posicao_celeste):
    """
    Calcula a distância radial de uma posição celeste.

    Parâmetros:
    posicao_celeste (array): Vetor de posição celeste.

    Retorna:
    float: Distância radial.
    """
    return np.linalg.norm(posicao_celeste)

def calcula_momento_angular(posicao_celeste, velocidade_celeste):
    """
    Calcula o momento angular.

    Parâmetros:
    posicao_celeste (array): Vetor de posição celeste.
    velocidade_celeste (array): Vetor de velocidade celeste.

    Retorna:
    array: Momento angular.
    """
    return np.cross(posicao_celeste, velocidade_celeste)

def calcula_excentricidade(velocidade_celeste, momento_angular, parametro_gravitacional, posicao_celeste, distancia_radial):
    """
    Calcula a excentricidade da órbita.

    Parâmetros:
    velocidade_celeste (array): Vetor de velocidade celeste.
    momento_angular (array): Vetor de momento angular.
    parametro_gravitacional (float): Parâmetro gravitacional.
    posicao_celeste (array): Vetor de posição celeste.
    distancia_radial (float): Distância radial.

    Retorna:
    array: Excentricidade da órbita.
    """
    return (np.cross(velocidade_celeste, momento_angular) / parametro_gravitacional - posicao_celeste / distancia_radial)

def calcula_parametros_orbitais(momento_angular, parametro_gravitacional, excentricidade):
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

    Parâmetros:
    posicao_celeste (array): Vetor de posição celeste.
    posicao_perifocal (array): Vetor de posição perifocal.
    distancia_radial (float): Distância radial.
    parametro (float): Parâmetro orbital.
    excentricidade (float): Excentricidade orbital.

    Retorna:
    float: Anomalia verdadeira.
    """
    cos_theta = (parametro - distancia_radial) / (excentricidade * distancia_radial)
    sin_theta = np.dot(posicao_celeste, posicao_perifocal) / (distancia_radial * parametro)
    return np.arctan2(sin_theta, cos_theta)

def calcula_tempo_periastro(epoca, excentricidade, anomalia_verdadeira, parametro_gravitacional, semi_eixo_maior):
    """
    Calcula o tempo de periastro.

    Parâmetros:
    epoca (float): Época inicial.
    excentricidade (float): Excentricidade orbital.
    anomalia_verdadeira (float): Anomalia verdadeira.
    parametro_gravitacional (float): Parâmetro gravitacional.
    semi_eixo_maior (float): Semi-eixo maior da órbita.

    Retorna:
    float: Tempo de periastro.
    """
    if excentricidade < 1:
        anomalia_excentrica = 2 * np.arctan(np.sqrt((1 - excentricidade) / (1 + excentricidade)) * np.tan(anomalia_verdadeira / 2))
        return epoca - (anomalia_excentrica - excentricidade * np.sin(anomalia_excentrica)) / np.sqrt(parametro_gravitacional / semi_eixo_maior ** 3)
    elif excentricidade == 1:
        parametro = semi_eixo_maior * (1 - excentricidade ** 2)  # Ajuste para órbita parabólica
        return epoca - (2 / 3) * ((parametro ** (3 / 2)) / np.sqrt(parametro_gravitacional)) * np.tan(anomalia_verdadeira / 2)
    else:
        anomalia_hiperbolica = 2 * np.arctanh(np.sqrt((excentricidade - 1) / (1 + excentricidade)) * np.tan(anomalia_verdadeira / 2))
        return epoca - (excentricidade * np.sinh(anomalia_hiperbolica) - anomalia_hiperbolica) / np.sqrt(-parametro_gravitacional / semi_eixo_maior ** 3)

def calcula_tempo_de_periastro_anomalia_verdadeira(momento_angular, excentricidade, distancia_radial, posicao_celeste, parametro, exc, epoca, parametro_gravitacional, semi_eixo_maior):
    """
    Calcula os elementos orbitais.

    Parâmetros:
    momento_angular (array): Vetor de momento angular.
    excentricidade (array): Vetor de excentricidade.
    distancia_radial (float): Distância radial.
    posicao_celeste (array): Vetor de posição celeste.
    parametro (float): Parâmetro orbital.
    exc (float): Excentricidade orbital.
    epoca (float): Época inicial.
    parametro_gravitacional (float): Parâmetro gravitacional.
    semi_eixo_maior (float): Semi-eixo maior da órbita.

    Retorna:
    tuple: Anomalia verdadeira e tempo de periastro.
    """
    posicao_perifocal = parametro * np.cross(momento_angular, excentricidade) / (np.linalg.norm(momento_angular) * exc)
    anomalia_verdadeira = calcula_anomalia_verdadeira(posicao_celeste, posicao_perifocal, distancia_radial, parametro, exc)
    tempo_periastro = calcula_tempo_periastro(epoca, exc, anomalia_verdadeira, parametro_gravitacional, semi_eixo_maior)
    return anomalia_verdadeira, tempo_periastro

def calcula_inclinacao_raan_argumento_de_periastro(momento_angular, excentricidade):
    """
    Calcula a inclinação, RAAN e argumento de periastro.

    Parâmetros:
    momento_angular (array): Vetor de momento angular.
    excentricidade (array): Vetor de excentricidade.

    Retorna:
    tuple: RAAN, inclinação e argumento de periastro.
    """
    momento_angular_normalizado = momento_angular / np.linalg.norm(momento_angular)
    vetor_referencia = np.array([0, 0, 1])
    linha_dos_nodos = np.cross(vetor_referencia, momento_angular_normalizado) / np.linalg.norm(np.cross(vetor_referencia, momento_angular_normalizado))
    raan = np.arctan2(linha_dos_nodos[1], linha_dos_nodos[0])
    inclinacao = np.arccos(np.dot(momento_angular_normalizado, vetor_referencia))
    excentricidade_normalizada = excentricidade / np.linalg.norm(excentricidade)
    cos_argumento_de_periastro = np.dot(excentricidade_normalizada, linha_dos_nodos)
    sin_argumento_de_periastro = np.dot(momento_angular_normalizado, np.cross(linha_dos_nodos, excentricidade_normalizada))
    argumento_de_periastro = np.arctan2(sin_argumento_de_periastro, cos_argumento_de_periastro)
    return raan, inclinacao, argumento_de_periastro

def determina_parametros_orbitais(epoca, parametro_gravitacional, posicao_celeste, velocidade_celeste):
    """
    Determina os parâmetros orbitais.

    Parâmetros:
    epoca (float): Época inicial.
    parametro_gravitacional (float): Parâmetro gravitacional.
    posicao_celeste (array): Vetor de posição celeste.
    velocidade_celeste (array): Vetor de velocidade celeste.

    Retorna:
    objeto: Órbita.
    """
    distancia_radial = calcula_distancia_radial(posicao_celeste)
    momento_angular = calcula_momento_angular(posicao_celeste, velocidade_celeste)
    excentricidade = calcula_excentricidade(velocidade_celeste, momento_angular, parametro_gravitacional, posicao_celeste, distancia_radial)
    exc = np.linalg.norm(excentricidade)
    h = np.linalg.norm(momento_angular)
    parametro, eixo = calcula_parametros_orbitais(h, parametro_gravitacional, exc)
    anomalia_verdadeira, tempo_periastro = calcula_tempo_de_periastro_anomalia_verdadeira(momento_angular, excentricidade, distancia_radial, posicao_celeste, parametro, exc, epoca, parametro_gravitacional, eixo)
    raan, inclinacao, argumento_de_periastro = calcula_inclinacao_raan_argumento_de_periastro(momento_angular, excentricidade)

    construtor_orbita = ConstrutorOrbita()
    orbita = (construtor_orbita
              .com_semi_eixo_maior(eixo)
              .com_excentricidade(exc)
              .com_inclinacao(inclinacao)
              .com_raan(raan)
              .com_arg_periastro(argumento_de_periastro)
              .com_anomalia_verdadeira(anomalia_verdadeira)
              .com_tempo_de_periastro(tempo_periastro)
              .construir())
    return orbita
