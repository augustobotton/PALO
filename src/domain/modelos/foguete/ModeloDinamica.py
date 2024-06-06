import numpy as np

from src.domain.modelos.manobras.parametros_manobra_adquire_gso import parametros_manobra_adquire_gso
from src.domain.modelos.planeta.ModeloGravidadePlanetaAxisimetrico import calcular_gravidade_axisimetric


def dinamica_foguete(vetor_tempo, vetor_de_estados_foguete, base_de_lancamento, planeta,
                     modelo_atmosferico, foguete):
	"""
	Função para a dinâmica de translação de um foguete com respeito ao referencial PCPF.
	Sistema de referência: aerodinâmico.
	Sistema de coordenadas: esférico.
	"""

	# Extrai os parâmetros do vetor de estado
	velocidade, azimute, phi, distancia_radial, latitude = vetor_de_estados_foguete[:5]

	if velocidade < 0:
		velocidade = 0.0001  # Evita velocidade negativa

	# Cálculo da massa e tração em função do tempo
	ft, massa, mu, epsl = foguete.modelo_propulsivo.propulsao_N_estagios(vetor_tempo,
	                                                                     vetor_de_estados_foguete)


	# Cálculo do modelo atmosférico
	altitude = distancia_radial - planeta.raio_equatorial
	T, _, _, densidade_do_ar, _, M, _, _, Kn, _, _, R = modelo_atmosferico.calcula(
		altitude, velocidade, planeta.tempo_longitude_celeste_nula, planeta.delta_temperatura_atm
	)

	# Cálculo do modelo aerodinâmico
	areas_de_referencia_para_calculo_do_arrasto, _, _ = (
		foguete.modelo_estrutural.calcula())
	tempos_de_separacao = foguete.modelo_propulsivo.tempos_de_separacao

	D, fy, L = foguete.modelo_aerodinamico.aerodinamica_multiplos_estagios(
		vetor_tempo, velocidade,areas_de_referencia_para_calculo_do_arrasto, tempos_de_separacao,
		densidade_do_ar)

	# Cálculo da gravidade
	gc, gd = calcular_gravidade_axisimetric(distancia_radial, latitude)

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

	# Saturação da altitude
	if altitude < 0:
		rp = 0
		deltap = 0
		lonp = 0
		Vp = 0
		Ap = 0
		phip = 0

	# Modelagem do trilho de lançamento
	altura_relativa = altitude - base_de_lancamento.altitude_base
	if altura_relativa <= base_de_lancamento.comprimento_trilho and vetor_tempo <= 10:
		Ap = 0
		phip = 0
	# TODO verificar onde e como chamar essa funcao provavelmente no foguete
	if not parametros_manobra_adquire_gso.achouapogeu:
		parametros_manobra_adquire_gso(vetor_tempo, massa, vetor_de_estados_foguete)

	# Derivada do vetor de estado
	Xp = [float(Vp), Ap, phip, rp, deltap, lonp]
	return Xp
