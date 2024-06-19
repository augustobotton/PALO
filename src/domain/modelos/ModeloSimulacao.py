import numpy as np
from scipy.integrate import solve_ivp

from src.domain.modelos.planeta.ModeloBaseDeLancamento import ConstrutorBaseDeLancamento
from src.domain.modelos.planeta.ModeloPlaneta import ConstrutorDePlanetas
class Simulacao:
	def __init__(self, planeta, base_de_lancamento, foguete, tempo_simulacao=150000,
	             velocidade_inicial=1, angulo_elevacao_inicial=80, inclinacao_orbita=5,
	             altitude_geo_sincrona=42.164140e6):
		self.tempo_simulacao = tempo_simulacao
		self.velocidade_inicial = velocidade_inicial
		self.angulo_elevacao_inicial = np.radians(angulo_elevacao_inicial)

		self.planeta = planeta
		self.base_de_lancamento = base_de_lancamento
		self.foguete = foguete

		self.inclinacao_orbita = np.radians(inclinacao_orbita)
		self.altitude_geo_sincrona = altitude_geo_sincrona

		self._configurar_parametros_planeta()
		self._configurar_parametros_base()
		self._configurar_parametros_orbita()
		self._calcular_condicoes_azimute()

	def _configurar_parametros_planeta(self):
		self.raio_equatorial = self.planeta.raio_equatorial
		self.velocidade_rotacao = self.planeta.velocidade_inercial_de_rotacao
		self.gravidade = self.planeta.gravidade_padrao_nivel_do_mar
		self.constante_gravitacional = self.planeta.mut
		self.j2 = self.planeta.J2
		self.j3 = self.planeta.J3
		self.j4 = self.planeta.J4

	def _configurar_parametros_base(self):
		self.altitude_inicial = self.base_de_lancamento.altitude_base_de_lancamento
		self.latitude_inicial = self.base_de_lancamento.latitude_base_de_lancamento
		self.longitude_inicial = self.base_de_lancamento.longitude_base_de_lancamento
		self.comprimento_trilho = self.base_de_lancamento.comprimento_do_trilho

	def _configurar_parametros_orbita(self):
		self.velocidade_geo_sincrona = np.sqrt(self.constante_gravitacional / self.altitude_geo_sincrona)

		self.massa_inicial = np.sum(self.foguete.mp) + np.sum(self.foguete.ms) + self.foguete.mL
		self.distancia_radial_inicial = self.raio_equatorial + self.altitude_inicial

	def simular(self):
		dinamica_foguete = self.foguete.dinamica_foguete
		condicoes_iniciais = [self.velocidade_inicial, self.azimute_inicial, self.angulo_elevacao_inicial,
		                      self.distancia_radial_inicial, self.latitude_inicial, self.longitude_inicial]
		opcoes_integracao = {'rtol': 1e-8, 'atol': 1e-10, 'max_step': 1}

		resposta_simulacao = solve_ivp(dinamica_foguete, (0, self.tempo_simulacao), y0=condicoes_iniciais,
		                               **opcoes_integracao)
		return resposta_simulacao.t, resposta_simulacao.y

	def _calcular_condicoes_azimute(self):
		y_t = np.cos(self.inclinacao_orbita) / np.cos(self.latitude_inicial)
		if abs(y_t) > 1:
			print(
				'Nao eh possivel atingir a inclinacao a partir da latitude inicial. Calculando a menor possivel')
			y_t = np.sign(y_t)

		self.azimute_final = np.arcsin(y_t)
		self.azimute_inicial = self._estimar_azimute_inicial()

		print(f'Condicao final de azimute de velocidade inercial (grau): {np.degrees(self.azimute_final)}')
		print(
			f'Condicao inicial de azimute de velocidade relativa (grau): {np.degrees(self.azimute_inicial)}')

	def _estimar_azimute_inicial(self):
		apogeu_transferencia = self.raio_equatorial + 250e3
		semi_eixo_maior_transferencia = (self.altitude_geo_sincrona + apogeu_transferencia) / 2
		velocidade_transferencia = np.sqrt(
			self.constante_gravitacional * (2 / apogeu_transferencia - 1 / semi_eixo_maior_transferencia))
		return np.arctan(np.tan(self.azimute_final) - (
					apogeu_transferencia * self.velocidade_rotacao * np.cos(self.latitude_inicial)) / (
					                 velocidade_transferencia * np.cos(self.azimute_final)))