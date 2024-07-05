class ModeloPlaneta:
    def __init__(self, delta_temperatura_atm: float, raio_equatorial: float,
                 velocidade_inercial_de_rotacao: float,
                 gravidade: float, mut: float, J2: float, J3: float, J4: float,
                 tempo_longitude_celeste_nula: float, modelo_atmosferico):
        self.delta_temperatura_atm = delta_temperatura_atm  # K - Delta T em relação à atmosfera padrão
        self.raio_equatorial = raio_equatorial  # m - Raio equatorial
        self.velocidade_inercial_de_rotacao = velocidade_inercial_de_rotacao  # Rad/s - Velocidade inercial de rotação
        self.gravidade = gravidade  # m/s^2 - Gravidade ao nível do mar
        self.mut = mut  # m^3/s^2 - Constante gravitacional
        self.J2 = J2  # Constante de Jeffery J2
        self.J3 = J3  # Constante de Jeffery J3
        self.J4 = J4  # Constante de Jeffery J4
        self.tempo_longitude_celeste_nula = tempo_longitude_celeste_nula  # s - Tempo de longitude celeste nula
        self.modelo_atmosferico = modelo_atmosferico

    def __str__(self):
        return f"ModeloPlaneta(ΔT_atm={self.delta_temperatura_atm}, R_eq={self.raio_equatorial}, " \
               f"ω={self.velocidade_inercial_de_rotacao}, g_0={self.gravidade}, " \
               f"μ={self.mut}, J2={self.J2}, J3={self.J3}, J4={self.J4}, t_λ=0={self.tempo_longitude_celeste_nula})"


