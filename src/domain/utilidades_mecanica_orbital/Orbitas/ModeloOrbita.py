import numpy as np

from src.domain.utilidades_mecanica_orbital.orbitalUtils import Converte


class Orbita:
	def __init__(self, semi_eixo_maior: float = None, excentricidade: float = None, inclinacao: float = None,
	             raan: float = None, arg_periastro: float = None, anomalia_verdadeira: float = None,
	             tempo_de_periastro: float = 0, parametro: float = None):
		"""
        Inicializa a classe Orbita com os parâmetros orbitais básicos.

        :param semi_eixo_maior: Semieixo maior da órbita (em metros).
        :param excentricidade: Excentricidade da órbita (adimensional).
        :param inclinacao: Inclinação da órbita.
        :param raan: Longitude do nó ascendente.
        :param arg_periastro: Argumento do periastro.
        :param anomalia_verdadeira: Anomalia verdadeira.
        """
		self.semi_eixo_maior = semi_eixo_maior
		self.excentricidade = excentricidade
		self.inclinacao = inclinacao
		self.raan = raan
		self.arg_periastro = arg_periastro
		self.anomalia_verdadeira = anomalia_verdadeira
		self.tempo_de_periastro = tempo_de_periastro
		# TODO criar metodo para overide mu
		self.parametro = parametro
		self.mu = 398600.4418e9  # Constante gravitacional da Terra, em m^3/s^2. É o padrão

	def calcular_parametro_orbital(self) -> float:
		"""
        Calcula o parâmetro orbital.

        :return: Parâmetro orbital.
        """
		return self.semi_eixo_maior * (1 - self.excentricidade ** 2)

	# def calcular_vetor_estado(self) -> tuple[np.ndarray, np.ndarray]:
	# 	"""
	#     Calcula o vetor de estado (posição e velocidade) a partir dos elementos orbitais.
	#
	#     :return: Vetores de posição e velocidade no sistema de coordenadas inerciais.
	#     """
	# 	p = self.calcular_parametro_orbital()
	# 	r = self._calcular_distancia_orbital(p)
	#
	# 	posicao_orbital = self._calcular_posicao_orbital(r)
	# 	velocidade_orbital = self._calcular_velocidade_orbital(p)
	#
	# 	matriz_rotacao = fdgfdgfdgfd
	#
	# 	posicao_inercial = np.dot(matriz_rotacao, posicao_orbital)
	# 	velocidade_inercial = np.dot(matriz_rotacao, velocidade_orbital)
	#
	# 	return posicao_inercial, velocidade_inercial

	def _calcular_distancia_orbital(self, p: float) -> float:
		"""
        Calcula a distância orbital.

        :param p: Parâmetro orbital.
        :return: Distância orbital.
        """
		return p / (1 + self.excentricidade * np.cos(self.anomalia_verdadeira))

	def _calcular_posicao_orbital(self, r: float) -> np.ndarray:
		"""
        Calcula a posição orbital.

        :param r: Distância orbital.
        :return: Vetor de posição no plano orbital.
        """
		x_orbital = r * np.cos(self.anomalia_verdadeira)
		y_orbital = r * np.sin(self.anomalia_verdadeira)
		z_orbital = 0
		return np.array([x_orbital, y_orbital, z_orbital])

	def _calcular_velocidade_orbital(self, p: float) -> np.ndarray:
		"""
        Calcula a velocidade orbital.

        :param p: Parâmetro orbital.
        :return: Vetor de velocidade no plano orbital.
        """
		h = np.sqrt(self.mu * p)
		vx_orbital = -self.mu / h * np.sin(self.anomalia_verdadeira)
		vy_orbital = self.mu / h * (self.excentricidade + np.cos(self.anomalia_verdadeira))
		vz_orbital = 0
		return np.array([vx_orbital, vy_orbital, vz_orbital])
