class ModeloPlaneta:
	def __init__(self, delta_temperatura_atm: float, raio_equatorial: float,
	             velocidade_inercial_de_rotacao: float,
	             gravidade_padrao_nivel_do_mar: float, mut: float, J2: float, J3: float, J4: float,
	             tempo_longitude_celeste_nula: float):
		self.delta_temperatura_atm = delta_temperatura_atm  # K - Delta T em relação à atmosfera padrão
		self.raio_equatorial = raio_equatorial  # m - Raio equatorial
		self.velocidade_inercial_de_rotacao = velocidade_inercial_de_rotacao  # Rad/s - Velocidade inercial de rotação
		self.gravidade_padrao_nivel_do_mar = gravidade_padrao_nivel_do_mar  # m/s^2 - Gravidade ao nível do mar
		self.mut = mut  # m^3/s^2 - Constante gravitacional
		self.J2 = J2  # Constante de Jeffery J2
		self.J3 = J3  # Constante de Jeffery J3
		self.J4 = J4  # Constante de Jeffery J4
		self.tempo_longitude_celeste_nula = tempo_longitude_celeste_nula  # s - Tempo de longitude celeste nula

	def __str__(self):
		return f"ModeloPlaneta(ΔT_atm={self.delta_temperatura_atm}, R_eq={self.raio_equatorial}, " \
		       f"ω={self.velocidade_inercial_de_rotacao}, g_0={self.gravidade_padrao_nivel_do_mar}, " \
		       f"μ={self.mut}, J2={self.J2}, J3={self.J3}, J4={self.J4}, t_λ=0={self.tempo_longitude_celeste_nula})"


class ConstrutorDePlanetas:
	def __init__(self):
		# Inicializa os atributos do construtor
		self.delta_temperatura_atm = None
		self.raio_equatorial = None
		self.velocidade_inercial_de_rotacao = None
		self.gravidade_padrao_nivel_do_mar = None
		self.mut = None
		self.J2 = None
		self.J3 = None
		self.J4 = None
		self.tempo_longitude_celeste_nula = None

	def com_delta_temperatura_atm(self, valor: float):
		self.delta_temperatura_atm = valor
		return self

	def com_raio_equatorial(self, valor: float):
		self.raio_equatorial = valor
		return self

	def com_velocidade_inercial_de_rotacao(self, valor: float):
		self.velocidade_inercial_de_rotacao = valor
		return self

	def com_gravidade_padrao_nivel_do_mar(self, valor: float):
		self.gravidade_padrao_nivel_do_mar = valor
		return self

	def com_mut(self, valor: float):
		self.mut = valor
		return self

	def com_J2(self, valor: float):
		self.J2 = valor
		return self

	def com_J3(self, valor: float):
		self.J3 = valor
		return self

	def com_J4(self, valor: float):
		self.J4 = valor
		return self

	def com_tempo_longitude_celeste_nula(self, valor: float):
		self.tempo_longitude_celeste_nula = valor
		return self

	def construir(self) -> ModeloPlaneta:
		# Valida os campos obrigatórios
		if None in (self.delta_temperatura_atm, self.raio_equatorial, self.velocidade_inercial_de_rotacao,
		            self.gravidade_padrao_nivel_do_mar, self.mut, self.J2, self.J3, self.J4,
		            self.tempo_longitude_celeste_nula):
			raise ValueError(
				"Todos os parâmetros devem ser fornecidos antes de construir o modelo do planeta.")

		# Retorna a instância de ModeloPlaneta configurada
		return ModeloPlaneta(
			self.delta_temperatura_atm,
			self.raio_equatorial,
			self.velocidade_inercial_de_rotacao,
			self.gravidade_padrao_nivel_do_mar,
			self.mut,
			self.J2,
			self.J3,
			self.J4,
			self.tempo_longitude_celeste_nula
		)


# Exemplo de uso
terra = (ConstrutorDePlanetas()
         .com_delta_temperatura_atm(10)
         .com_raio_equatorial(6378.1370e3)
         .com_velocidade_inercial_de_rotacao(7.2921150e-5)
         .com_gravidade_padrao_nivel_do_mar(9.80665)
         .com_mut(3.986004418e14)
         .com_J2(0.00108263)
         .com_J3(-0.00000254)
         .com_J4(-0.00000161)
         .com_tempo_longitude_celeste_nula(0)
         .construir())

print(terra)
