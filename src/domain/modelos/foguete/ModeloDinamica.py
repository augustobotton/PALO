import numpy as np

from src.domain.modelos.planeta.ModeloGravidadePlanetaAxisimetrico import calcular_gravidade_axisimetrico


def dinamica_foguete(vetor_tempo, vetor_de_estados_foguete, base_de_lancamento, planeta,
                     foguete, parametros_apogeu, orbita_alvo):
    """
	Função para a dinâmica de translação de um foguete com respeito ao referencial PCPF.
	Sistema de referência: aerodinâmico.
	Sistema de coordenadas: esférico.
	"""

    velocidade, azimute, phi, distancia_radial, latitude, longitude = vetor_de_estados_foguete[:6]

    if velocidade < 0:
        velocidade = 0.0001  # Evita velocidade negativa

    # Cálculo da massa e tração em função do tempo
    ft, massa, mu, epsl = foguete.modelo_propulsivo.propulsao_n_estagios(vetor_tempo,
                                                                         vetor_de_estados_foguete,
                                                                         foguete.modelo_estrutural)
    # Cálculo do modelo atmosférico
    altitude = distancia_radial - planeta.raio_equatorial

    # Cálculo do modelo aerodinâmico
    T, _, _, densidade_do_ar, _, M, _, _, Kn, _, _, R = planeta.modelo_atmosferico.calcula(
        altitude, velocidade, foguete.modelo_estrutural.comprimento_caracteristico, planeta.delta_temperatura_atm
    )

    foguete.modelo_aerodinamico.atualizar_parametros(altitude=altitude, numero_de_knudsen=Kn, numero_de_mach=M,
                                                     temperatura=T, constante_do_gas_ideal=R, velocidade=velocidade)
    tempos_de_separacao = foguete.modelo_propulsivo.tempos_de_separacao
    # Cálculo da força de arrasto
    D, fy, L = foguete.modelo_aerodinamico.aerodinamica_multiplos_estagios(
        vetor_tempo, velocidade, foguete.modelo_estrutural.areas_de_referencia_para_calculo_do_arrasto,
        tempos_de_separacao,
        densidade_do_ar)

    # Cálculo da gravidade
    gc, gd = calcular_gravidade_axisimetrico(distancia_radial, latitude, planeta)

    # Equações de cinemática de translação
    rp = velocidade * np.sin(phi)

    deltap = (velocidade / distancia_radial) * np.cos(phi) * np.cos(azimute)

    lonp = (velocidade * np.cos(phi) * np.sin(azimute)) / (distancia_radial * np.cos(latitude))

    # Equações de dinâmica de translação
    Vp = (1 / massa) * (
            ft * np.cos(epsl) * np.cos(mu) - D - massa * gc * np.sin(phi) +
            massa * gd * np.cos(phi) * np.cos(azimute) -
            massa * planeta.velocidade_inercial_de_rotacao ** 2 * distancia_radial * np.cos(latitude) * (
                    np.cos(phi) * np.cos(azimute) * np.sin(latitude) - np.sin(phi) * np.cos(latitude))
    )

    Ap = (1 / (massa * velocidade * np.cos(phi))) * (
            massa * (velocidade ** 2 / distancia_radial) * np.cos(phi) ** 2 * np.sin(azimute) * np.tan(
        latitude) +
            ft * np.sin(mu) + fy - massa * gd * np.sin(azimute) +
            massa * planeta.velocidade_inercial_de_rotacao ** 2 * distancia_radial * np.sin(azimute) *
            np.sin(
                latitude) * np.cos(latitude) -
            2 * massa * planeta.velocidade_inercial_de_rotacao * velocidade * (
                    np.sin(phi) * np.cos(azimute) * np.cos(latitude) - np.cos(phi) * np.sin(latitude))
    )

    phip = (1 / (massa * velocidade)) * (
            massa * (velocidade ** 2 / distancia_radial) * np.cos(phi) + ft * np.sin(epsl) * np.cos(mu) + L -
            massa * gc * np.cos(phi) - massa * gd * np.sin(phi) * np.cos(azimute) +
            massa * planeta.velocidade_inercial_de_rotacao ** 2 * distancia_radial * np.cos(latitude) * (
                    np.sin(phi) * np.cos(azimute) * np.sin(latitude) + np.cos(phi) * np.cos(latitude)) +
            2 * massa * planeta.velocidade_inercial_de_rotacao * velocidade * np.sin(azimute) * np.cos(
        latitude)
    )

    # Saturação da altitude - caso negativo o foguete fica na superfície
    if altitude < 0:
        rp = 0
        deltap = 0
        lonp = 0
        Vp = 0
        Ap = 0
        phip = 0

    # Modelagem do trilho de lançamento
    altura_relativa = altitude - base_de_lancamento.altitude_base_de_lancamento
    if altura_relativa <= base_de_lancamento.comprimento_do_trilho and vetor_tempo <= 10:
        Ap = 0
        phip = 0

    parametros_apogeu.parametros_manobra_adquire_gso(vetor_tempo, massa, vetor_de_estados_foguete, orbita_alvo,
                                                     foguete.modelo_propulsivo,
                                                     planeta)

    # Derivada do vetor de estado
    Xp = np.array([float(Vp), Ap, phip, rp, deltap, lonp])
    return Xp
