class ParametrosPlanetarios:
    def __init__(self, delta_temperatura_atm, raio_equatorial, velocidade_inercial_de_rotacao,
                 gravidade_padrao_nivel_do_mar, mut, J2, J3, J4, tempo_longitude_celeste_nula):
        self.delta_temperatura_atm = delta_temperatura_atm  # K - Delta T em relação à atmosfera padrão
        self.raio_equatorial = raio_equatorial  # m - Raio equatorial
        self.velocidade_inercial_de_rotacao = velocidade_inercial_de_rotacao
        self.gravidade_padrao_nivel_do_mar = gravidade_padrao_nivel_do_mar
        self.mut = mut  # m^3/s^-2
        self.J2 = J2  # Constante de Jeffery J2
        self.J3 = J3  # Constante de Jeffery J3
        self.J4 = J4  # Constante de Jeffery J4
        self.tempo_longitude_celeste_nula = tempo_longitude_celeste_nula


class ParametrosPlanetariosBuilder:
    def __init__(self):
        self.delta_temperatura_atm = None
        self.raio_equatorial = None
        self.velocidade_inercial_de_rotacao = None
        self.gravidade_padrao_nivel_do_mar = None
        self.mut = None
        self.J2 = None
        self.J3 = None
        self.J4 = None
        self.tempo_longitude_celeste_nula = None

    def com_delta_temperatura_atm(self, valor):
        self.delta_temperatura_atm = valor
        return self

    def com_raio_equatorial(self, valor):
        self.raio_equatorial = valor
        return self

    def com_velocidade_inercial_de_rotacao(self, valor):
        self.velocidade_inercial_de_rotacao = valor
        return self

    def com_gravidade_padrao_nivel_do_mar(self, valor):
        self.gravidade_padrao_nivel_do_mar = valor
        return self

    def com_mut(self, valor):
        self.mut = valor
        return self

    def com_J2(self, valor):
        self.J2 = valor
        return self

    def com_J3(self, valor):
        self.J3 = valor
        return self

    def com_J4(self, valor):
        self.J4 = valor
        return self

    def com_tempo_longitude_celeste_nula(self, valor):
        self.tempo_longitude_celeste_nula = valor
        return self

    def construir(self):
        return ParametrosPlanetarios(
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


# Exemplo de uso do Builder para criar um objeto ParametrosPlanetarios
parametros = ParametrosPlanetariosBuilder() \
    .com_delta_temperatura_atm(10) \
    .com_raio_equatorial(6378.1370e3) \
    .com_velocidade_inercial_de_rotacao(7.2921150e-5) \
    .com_gravidade_padrao_nivel_do_mar(9.80665) \
    .com_mut(3.986004418e14) \
    .com_J2(0.00108263) \
    .com_J3(-0.00000254) \
    .com_J4(-0.00000161) \
    .com_tempo_longitude_celeste_nula(0) \
    .construir()
