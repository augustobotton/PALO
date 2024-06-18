import numpy as np

from src.domain.utilidades_mecanica_orbital.Orbitas.ModeloOrbita import Orbita


class Manobras:
	def __init__(self, orbita):
		"""
		Inicializa a classe Manobras com uma instância da classe Orbita.

		:param orbita: Instância da classe Orbita.
		"""
		self.orbita = orbita

	def delta_v_mudanca_de_velocidade(self, delta_v):
		"""
		Aplica uma mudança de velocidade (delta-v) à órbita.

		:param delta_v: Vetor delta-v (em m/s).
		"""
		posicao, velocidade = self.orbita.vetor_estado()
		nova_velocidade = velocidade + delta_v
		self._atualizar_orbita(posicao, nova_velocidade)

	def mudanca_de_altitude(self, nova_altitude):
		"""
		Altera a altitude da órbita (supõe órbita circular para simplicidade).

		:param nova_altitude: Nova altitude desejada (em metros).
		"""
		mu = 398600.4418e9  # Constante gravitacional da Terra, em m^3/s^2
		raio_terra = 6371e3  # Raio da Terra, em metros

		nova_distancia = raio_terra + nova_altitude
		nova_velocidade_orbital = np.sqrt(mu / nova_distancia)

		posicao, _ = self.orbita.vetor_estado()
		nova_velocidade = np.array([0, nova_velocidade_orbital, 0])

		self._atualizar_orbita(posicao, nova_velocidade)

	def _atualizar_orbita(self, posicao, velocidade):
		"""
		Atualiza os elementos orbitais com base nos novos vetores de posição e velocidade.

		:param posicao: Vetor de posição (em metros).
		:param velocidade: Vetor de velocidade (em m/s).
		"""
		# Cálculo dos novos parâmetros orbitais a partir dos vetores de posição e velocidade.
		mu = 398600.4418e9  # Constante gravitacional da Terra, em m^3/s^2
		h_vector = np.cross(posicao, velocidade)
		h = np.linalg.norm(h_vector)

		r = np.linalg.norm(posicao)
		v = np.linalg.norm(velocidade)

		semi_eixo_maior = 1 / (2 / r - v ** 2 / mu)
		excentricidade_vector = (1 / mu) * (np.cross(velocidade, h_vector) - mu * posicao / r)
		excentricidade = np.linalg.norm(excentricidade_vector)

		inclinacao = np.arccos(h_vector[2] / h)
		raan = np.arctan2(h_vector[0], -h_vector[1])

		nodo_ascendente_vector = np.cross([0, 0, 1], h_vector)
		nodo_ascendente = nodo_ascendente_vector / np.linalg.norm(nodo_ascendente_vector)

		arg_perigeu = np.arccos(np.dot(nodo_ascendente, excentricidade_vector) / excentricidade)
		if excentricidade_vector[2] < 0:
			arg_perigeu = 2 * np.pi - arg_perigeu

		anomalia_verdadeira = np.arccos(np.dot(excentricidade_vector, posicao) / (excentricidade * r))
		if np.dot(posicao, velocidade) < 0:
			anomalia_verdadeira = 2 * np.pi - anomalia_verdadeira

		self.orbita.semi_eixo_maior = semi_eixo_maior
		self.orbita.excentricidade = excentricidade
		self.orbita.inclinacao = inclinacao
		self.orbita.raan = raan
		self.orbita.arg_periastro = arg_perigeu
		self.orbita.anomalia_verdadeira = anomalia_verdadeira


# Exemplo de uso
orbita_inicial = Orbita(semi_eixo_maior=7000e3, excentricidade=0.001, inclinacao=98.7, raan=120,
                        arg_periastro=45, anomalia_verdadeira=10)
manobras = Manobras(orbita_inicial)

# Aplicar uma mudança de velocidade
delta_v = np.array([0, 0.1, 0])  # Mudança de velocidade de 0.1 m/s na direção y
manobras.delta_v_mudanca_de_velocidade(delta_v)

# Alterar a altitude da órbita
nova_altitude = 800e3  # Nova altitude de 800 km
manobras.mudanca_de_altitude(nova_altitude)

# Verificar os novos parâmetros orbitais
print("Novo Semieixo Maior:", orbita_inicial.semi_eixo_maior)
print("Nova Excentricidade:", orbita_inicial.excentricidade)
print("Nova Inclinação:", np.degrees(orbita_inicial.inclinacao))
print("Novo RAAN:", np.degrees(orbita_inicial.raan))
print("Novo Argumento do Perigeu:", np.degrees(orbita_inicial.arg_periastro))
print("Nova Anomalia Verdadeira:", np.degrees(orbita_inicial.anomalia_verdadeira))