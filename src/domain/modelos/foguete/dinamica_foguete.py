import numpy as np

from src.domain.manobras.parametros_manobra_adquire_gso import parametros_manobra_adquire_gso
from src.domain.modelos.planeta.ModeloAtmosferico import ModeloAtmosferico
from src.domain.modelos.planeta.ModeloBaseDeLancamento import ConstrutorBaseDeLancamento
from src.domain.modelos.planeta.ModeloGravidadePlanetaAxisimetrico import grav_axisimetrico
from src.domain.modelos.planeta.ModeloPlaneta import ConstrutorDePlanetas


def dinamica_foguete(vetor_tempo, vetor_de_estados_foguete, base_de_lancamento):
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
    planeta = ConstrutorDePlanetas() \
        .com_delta_temperatura_atm(10) \
        .com_raio_equatorial(6378.1370e3) \
        .com_velocidade_inercial_de_rotacao(7.2921150e-5) \
        .com_gravidade_padrao_nivel_do_mar(9.80665) \
        .com_mut(3.986004418e14) \
        .com_J2(0.00108263) \
        .com_J3(-0.00000254) \
        .com_J4(-0.00000161) \
        .com_tempo_longitude_celeste_nula(0) \
        .construir()

    ## Entrada de constantes por variáveis globais
    we = planeta.velocidade_inercial_de_rotação_da_terra
    Re = planeta.raio_equatorial
    lc = planeta.tempo_longitude_celeste_nula
    dT = planeta.delta_temperatura_atm

    base_de_lancamento = (ConstrutorBaseDeLancamento().com_altitude_base(0)
                          .com_latitude_inicial(None)
                          .com_longitude_inicial(None)
                          .com_comprimento_trilho(10)
                          .construir())
    achouApogeu = base_de_lancamento.achouApogeu

    ## Vetor de estado
    velocidade = vetor_de_estados_foguete[0]
    azimute_vetor_velocidade_relativa = vetor_de_estados_foguete[1]
    phi = vetor_de_estados_foguete[2]
    distancia_radial_centro_da_terra = vetor_de_estados_foguete[3]
    latitude_em_relacao_ao_plano_equatorial = vetor_de_estados_foguete[4]
    # lon = X[5]

    if velocidade < 0:
        velocidade = 0.0001  # Evita velocidade negativa

    ## Função para cálculo da massa e tração em função do tempo
    # Depende da taxa de queima de propelente, tempo de queima de propelente,
    # impulso específico e separação dos estágios
    ft, m, mu, epsl = propulsao_N_estagios(vetor_tempo, vetor_de_estados_foguete)

    ## Função para cálculo do modelo atmosférico
    h = distancia_radial_centro_da_terra - Re  # Altitude
    modelo_atmosferico = ModeloAtmosferico()
    T, _, _, rho, _, M, _, _, Kn, _, _, R = modelo_atmosferico.calcula(h, velocidade, lc, dT)

    ## Função para cálculo do modelo aerodinâmico
    # Depende da altitude e velocidade
    D, fy, L = aerodinamica_multiplos_estagios(vetor_tempo, velocidade, h, M, Kn, T, rho, R)

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
    H = h - base_de_lancamento.altitude_base_de_lancamento  # Altura
    if ((H <= base_de_lancamento.comprimento_do_trilho) and (
            vetor_tempo <= 10)):  # Verifica se a altura eh menor que l_trilho nos primeiros segundos da simulacao
        Ap = 0
        phip = 0  # Anula as derivadas dos angulos de orientacao da velocidade

    if achouApogeu == 0:
        parametros_manobra_adquire_gso(vetor_tempo, m, vetor_de_estados_foguete)
    ## Derivada do vetor de estado
    Xp = [float(Vp), Ap, phip, rp, deltap, lonp]
    return Xp
