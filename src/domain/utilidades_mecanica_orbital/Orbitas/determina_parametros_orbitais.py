import numpy as np

from src.domain.utilidades_mecanica_orbital.Orbitas.ConstrutorOrbita import ConstrutorOrbita


def determina_parametros_orbitais(t0, mu, posicao_celeste, velocidade_celeste):
	# Distância radial ao primário no instante observado
	distancia_radial_observada = np.linalg.norm(posicao_celeste)

	# Vetor quantidade de movimento angular específica no referencial celeste
	hc = np.cross(posicao_celeste, velocidade_celeste)

	excentricidade_celeste = (np.cross(velocidade_celeste, hc) / mu - posicao_celeste /
	                          distancia_radial_observada)

	e = np.linalg.norm(excentricidade_celeste)

	h = np.linalg.norm(hc)

	# Semi-eixo maior e parâmetro da órbita

	p = h ** 2 / mu
	a = p / (1 - e ** 2)

	# Vetor parâmetro no referencial celeste
	pc = p * np.cross(hc, excentricidade_celeste) / (h * e)

	# Anomalia Verdadeira
	cos_theta = (p - distancia_radial_observada) / (e * distancia_radial_observada)
	sin_theta = np.dot(posicao_celeste, pc) / (distancia_radial_observada * p)
	anomalia_verdadeira = np.arctan2(sin_theta, cos_theta)


	if e < 1:
		# Órbita elíptica
		E = 2 * np.arctan(np.sqrt((1 - e) / (1 + e)) * np.tan(anomalia_verdadeira / 2))
		tempo_periastro = t0 - (E - e * np.sin(E)) / np.sqrt(mu / a ** 3)
	elif e == 1:
		# Órbita parabólica
		tempo_periastro = -((np.tan(anomalia_verdadeira / 2)) ** 3 + 3 * np.tan(anomalia_verdadeira / 2)) / (
				mu / p ** 3) ** (1 / 6)
	else:
		# Órbita hiperbólica
		H = 2 * np.arctanh(np.sqrt((e - 1) / (1 + e)) * np.tan(anomalia_verdadeira / 2))
		tempo_periastro = t0 - (e * np.sinh(H) - H) / np.sqrt(-mu / a ** 3)

	# Vetor unitário ao longo do vetor h (no sistema celeste)
	ih = hc / h
	#  Vetor unitário ao longo da linha dos nodos (no sistema celeste)
	Kc = np.array([0, 0, 1])
	nc = np.cross(Kc, ih) / np.linalg.norm(np.cross(Kc, ih))

	raan = np.arctan2(nc[1], nc[0])
	i = np.arccos(np.dot(ih, Kc))
	ie = excentricidade_celeste / e

	cos_omega = np.dot(ie, nc)
	sin_omega = np.dot(ih, np.cross(nc, ie))
	argumento_periastro = np.arctan2(sin_omega, cos_omega)
	criar_orbita = ConstrutorOrbita()
	orbita = criar_orbita.com_semi_eixo_maior(a).com_excentricidade(e).com_inclinacao(i).com_raan(
		raan).com_arg_periastro(argumento_periastro).com_anomalia_verdadeira(
		anomalia_verdadeira).com_tempo_de_periastro(tempo_periastro).construir()

	return orbita

