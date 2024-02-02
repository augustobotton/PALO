import numpy as np
from aerodinamica_N_estagios import aerodinamica_multiplos_estagios
from parametros_manobra_adquire_gso import parametros_manobra_adquire_gso

from propulsao_N_estagios import propulsao_N_estagios
from src.domain.modelos.planeta.ModeloAtmosferico import ModeloAtmosferico
from src.domain.modelos.planeta.grav_axisimetrico import grav_axisimetrico


def dinamica_foguete(t, X, terra, base_de_lancamento):
    """
    Função para a dinâmica de translação de um foguete com respeito ao referencial PCPF
    Sistema de referência: aerodinâmico
    Sistema de coordenadas: esférico
    """
    # Entradas:
    # t: tempo (s)
    # X: vetor de estado
    # V = X(1) (m/s): módulo do vetor velocidade relativa com respeito ao planeta girante
    # A = X(2) (rad): ângulo de azimute do vetor velocidade relativa com respeito ao eixo z (que aponta para o norte) do sistema uen.
    # phi = X(3) (rad): ângulo de elevação do vetor velocidade relativa com respeito ao horizonte local (plano yz do referencial uen)
    # r = X(4) (m): distância radial até o centro do planeta
    # delta = X(5) (rad): latitude com respeito ao plano equatorial do planeta
    # lon = X(6) (rad): longitude planetaria
    # Saída:
    # Xp: derivada do vetor de estado X

    ## Entrada de constantes por variáveis globais
    we = terra.velocidade_inercial_de_rotação_da_terra
    Re = terra.raio_equatorial
    lc = terra.tempo_longitude_celeste_nula
    dT = terra.delta_temperatura_atm
    h0 = parametros.altitude_base_de_lancamento
    l_trilho = parametros.comprimento_do_trilho
    achouApogeu = parametros.achouApogeu

    ## Vetor de estado
    velocidade = X[0]
    azimute_vetor_velocidade_relativa = X[1]
    phi = X[2]
    distancia_radial_centro_da_terra = X[3]
    latitude_em_relacao_ao_plano_equatorial = X[4]
    # lon = X[5]

    if velocidade < 0:
        velocidade = 0.0001;  # Evita velocidade negativa

    ## Função para cálculo da massa e tração em função do tempo
    # Depende da taxa de queima de propelente, tempo de queima de propelente,
    # impulso específico e separação dos estágios
    ft, m, mu, epsl = propulsao_N_estagios(t, X)

    ## Função para cálculo do modelo atmosférico
    h = distancia_radial_centro_da_terra - Re  # Altitude
    modelo_atmosferico = ModeloAtmosferico()
    T, _, _, rho, _, M, _, _, Kn, _, _, R = modelo_atmosferico.calcula(h, velocidade, lc, dT)

    ## Função para cálculo do modelo aerodinâmico
    # Depende da altitude e velocidade
    D, fy, L = aerodinamica_multiplos_estagios(t, velocidade, h, M, Kn, T, rho, R)

    ## Calculo da gravidade
    # Função para cálculo do modelo gravitacional
    gc, gd = grav_axisimetrico(distancia_radial_centro_da_terra, latitude_em_relacao_ao_plano_equatorial)

    ## Equações de cinemática de translação
    rp = velocidade * np.sin(phi)

    deltap = (velocidade / distancia_radial_centro_da_terra) * np.cos(phi) * np.cos(azimute_vetor_velocidade_relativa)

    lonp = (velocidade * np.cos(phi) * np.sin(azimute_vetor_velocidade_relativa)) / (
            distancia_radial_centro_da_terra * np.cos(latitude_em_relacao_ao_plano_equatorial))

    ## Equações de dinâmica de translação
    Vp = (1 / m) * (ft * np.cos(epsl) * np.cos(mu) - D - m * gc * np.sin(phi) + m * gd * np.cos(phi) * np.cos(
        azimute_vetor_velocidade_relativa) - \
                    m * we ** 2 * distancia_radial_centro_da_terra * np.cos(latitude_em_relacao_ao_plano_equatorial) * (
                            np.cos(phi) * np.cos(azimute_vetor_velocidade_relativa) * np.sin(
                        latitude_em_relacao_ao_plano_equatorial) - np.sin(phi) * np.cos(
                        latitude_em_relacao_ao_plano_equatorial)))
    # Vp = (1 / m) * (ft * np.cos(epsl) * np.cos(mu) - D - m * gc * np.sin(phi) + m * gd * np.cos(phi) * np.cos(A) - m * we**2 * r * np.cos(delta) * (np.cos(phi) * np.cos(A) * np.sin(delta) - np.sin(phi) * np.cos(delta)))
    Ap = (1 / (m * velocidade * (np.cos(phi)))) * (
            m * (velocidade ** 2 / distancia_radial_centro_da_terra) * np.cos(phi) ** 2 * np.sin(
        azimute_vetor_velocidade_relativa) * np.tan(latitude_em_relacao_ao_plano_equatorial) + ft * np.sin(
        mu) + fy - m * gd * np.sin(azimute_vetor_velocidade_relativa) + \
            m * we ** 2 * distancia_radial_centro_da_terra * np.sin(azimute_vetor_velocidade_relativa) * np.sin(
        latitude_em_relacao_ao_plano_equatorial) * np.cos(
        latitude_em_relacao_ao_plano_equatorial) - 2 * m * we * velocidade * (
                    np.sin(phi) * np.cos(azimute_vetor_velocidade_relativa) * np.cos(
                latitude_em_relacao_ao_plano_equatorial) - np.cos(phi) * np.sin(
                latitude_em_relacao_ao_plano_equatorial)));
    # Ap = (1 / (m * V * np.cos(phi))) * (m * (V**2 / r) * np.cos(phi)**2 * np.sin(A) * np.tan(delta) + ft * np.sin(mu) + fy - m * gd * np.sin(A) + m * we**2 * r * np.sin(A) * np.sin(delta) * np.cos(delta) - 2 * m * we * V * (np.sin(phi) * np.cos(A) * np.cos(delta) - np.cos(phi) * np.sin(delta)))
    phip = 1 * (1 / (m * velocidade)) * (
            m * (velocidade ** 2 / distancia_radial_centro_da_terra) * np.cos(phi) + ft * np.sin(epsl) * np.cos(
        mu) + L - m * gc * np.cos(phi) - m * gd * np.sin(phi) * np.cos(azimute_vetor_velocidade_relativa) \
            + m * we ** 2 * distancia_radial_centro_da_terra * np.cos(latitude_em_relacao_ao_plano_equatorial) * (
                    np.sin(phi) * np.cos(azimute_vetor_velocidade_relativa) * np.sin(
                latitude_em_relacao_ao_plano_equatorial) + np.cos(phi) * np.cos(
                latitude_em_relacao_ao_plano_equatorial)) + 2 * m * we * velocidade * np.sin(
        azimute_vetor_velocidade_relativa) * np.cos(latitude_em_relacao_ao_plano_equatorial));
    # phip = (1 / (m * V)) * (m * (V**2 / r) * np.cos(phi) + ft * np.sin(epsl) * np.cos(mu) + L - m * gc * np.cos(phi) - m * gd * np.sin(phi) * np.cos(A) + m * we**2 * r * np.cos(delta) * (np.sin(phi) * np.cos(A) * np.sin(delta) + np.cos(phi) * np.cos(delta)) + 2 * m * we * V * np.sin(A) * np.cos(delta))

    ## Saturação da altitude
    if h < 0:  # Altitude negativa não é permitida
        # Mantém as derivadas nulas
        rp = 0
        deltap = 0
        lonp = 0
        Vp = 0
        Ap = 0
        phip = 0

    ## Modela o trilho de lancamento
    H = h - h0  # Altura
    if ((H <= l_trilho) and (
            t <= 10)):  # Verifica se a altura eh menor que l_trilho nos primeiros segundos da simulacao
        Ap = 0
        phip = 0  # Anula as derivadas dos angulos de orientacao da velocidade

    if achouApogeu == 0:
        parametros_manobra_adquire_gso(t, m, X)
    ## Derivada do vetor de estado
    Xp = [float(Vp), Ap, phip, rp, deltap, lonp]
    return Xp
