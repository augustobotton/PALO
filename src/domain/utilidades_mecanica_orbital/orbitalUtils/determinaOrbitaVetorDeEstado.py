import numpy as np

from src.domain.utilidades_mecanica_orbital.Orbitas.ModeloOrbita import Orbita


class UtilidadesOrbitais:
	@staticmethod
	def matriz_rotacao_orbital_inercial(raan: float, arg_perigeu: float, inclinacao: float) -> np.ndarray:
		"""
		Cria a matriz de rotação do plano orbital para o sistema de coordenadas inerciais.

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

	@staticmethod
	def det_orbita(t0: float, rc0: np.ndarray, vc0: np.ndarray, mu: float) -> 'Orbita':
		"""
		Determina os parâmetros orbitais a partir de uma observação de posição e velocidade.

		:param t0: (s) Tempo em que a observação foi feita.
		:param rc0: (m ou km) Vetor posição relativa de m2 com respeito a m1 escrito no referencial celeste.
		:param vc0: (m/s ou km/s) Vetor velocidade relativa escrito no referencial celeste.
		:param mu: (m^3/s^2 ou km^3/s^2) Parâmetro gravitacional padrão do corpo m1.
		:return: Objeto da classe Orbita com os parâmetros orbitais calculados.
		"""
		r0 = np.linalg.norm(rc0)
		hc = np.cross(rc0, vc0)
		ec = np.cross(vc0, hc) / mu - rc0 / r0
		e = np.linalg.norm(ec)
		h = np.linalg.norm(hc)
		p = h ** 2 / mu
		a = p / (1 - e ** 2)
		pc = p * np.cross(hc, ec) / (h * e)
		costheta = (p - r0) / (e * r0)
		sintheta = np.dot(rc0, pc) / (r0 * p)
		theta0 = np.arctan2(sintheta, costheta)

		if 0 <= e < 1:
			tipo = 'e'
		elif e == 1:
			tipo = 'p'
		else:
			tipo = 'h'

		if tipo == 'e':
			n = np.sqrt(mu / a ** 3)
			E0 = 2 * np.arctan(np.sqrt((1 - e) / (1 + e)) * np.tan(theta0 / 2))
			tau = t0 - (E0 - e * np.sin(E0)) / n
		elif tipo == 'p':
			tau = t0 - (np.tan(theta0 / 2) ** 3 / 6 + np.tan(theta0 / 2) / 2) / np.sqrt(mu / p ** 3)
		else:
			H0 = 2 * np.arctanh(np.sqrt((e - 1) / (1 + e)) * np.tan(theta0 / 2))
			tau = t0 - (e * np.sinh(H0) - H0) / ((e ** 2 - 1) ** (3 / 2) * np.sqrt(mu / p ** 3))

		Kc = np.array([0, 0, 1])
		nc = np.cross(Kc, hc) / np.linalg.norm(np.cross(Kc, hc))
		OMEGA = np.arctan2(nc[1], nc[0])

		if abs(OMEGA) > 1e-3:
			i = np.arctan2(hc[0] / np.sin(OMEGA), hc[2])
		else:
			i = np.arctan2(-hc[1] / np.cos(OMEGA), hc[2])

		ie = ec / e
		ih = hc / h
		cosomega = np.dot(ie, nc)
		sinomega = np.dot(ih, np.cross(nc, ie))
		omega = np.arctan2(sinomega, cosomega)

		# Criar objeto Orbita com os parâmetros calculados
		return Orbita(a, e, np.degrees(i), np.degrees(OMEGA), np.degrees(omega), np.degrees(theta0))