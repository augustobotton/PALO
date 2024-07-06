class ModeloPlaneta:
    """
    Representa as propriedades físicas e atmosféricas de um planeta para simulações de dinâmica de voo.

    Attributes:
        delta_temperatura_atm (float): Delta de temperatura em relação à atmosfera padrão em Kelvin.
        raio_equatorial (float): Raio equatorial do planeta em metros.
        velocidade_inercial_de_rotacao (float): Velocidade inercial de rotação do planeta em radianos por segundo.
        gravidade (float): Aceleração devido à gravidade ao nível do mar em m/s^2.
        mut (float): Constante gravitacional do planeta em m^3/s^2.
        J2 (float): Coeficiente do termo J2 no desenvolvimento em série do potencial gravitacional.
        J3 (float): Coeficiente do termo J3.
        J4 (float): Coeficiente do termo J4.
        tempo_longitude_celeste_nula (float): Tempo de longitude celeste nula em segundos.
        modelo_atmosferico (object): Modelo atmosférico do planeta.
    """

    def __init__(self, delta_temperatura_atm: float, raio_equatorial: float,
                 velocidade_inercial_de_rotacao: float, gravidade: float, mut: float,
                 J2: float, J3: float, J4: float, tempo_longitude_celeste_nula: float,
                 modelo_atmosferico, massa: float):
        self.delta_temperatura_atm = delta_temperatura_atm
        self.raio_equatorial = raio_equatorial
        self.velocidade_inercial_de_rotacao = velocidade_inercial_de_rotacao
        self.gravidade = gravidade
        self.mut = mut
        self.J2 = J2
        self.J3 = J3
        self.J4 = J4
        self.tempo_longitude_celeste_nula = tempo_longitude_celeste_nula
        self.modelo_atmosferico = modelo_atmosferico
        self.massa = massa

    def __str__(self):
        return (f"ModeloPlaneta(ΔT_atm={self.delta_temperatura_atm} K, R_eq={self.raio_equatorial} m, "
                f"ω={self.velocidade_inercial_de_rotacao} rad/s, g_0={self.gravidade} m/s², "
                f"μ={self.mut} m³/s², J2={self.J2}, J3={self.J3}, J4={self.J4}, "
                f"t_λ=0={self.tempo_longitude_celeste_nula} s)")
