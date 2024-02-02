import numpy as np


class BaseDeLancamento:
    def __init__(self):
        self.altitude_base_de_lancamento = None
        self.latitude_base_de_lancamento = None
        self.longitude_base_de_lancamento = None
        self.comprimento_do_trilho = None

    def __str__(self):
        return f"Altitude da base de lançamento: {self.altitude_base_de_lancamento} m, " \
               f"Latitude inicial: {self.latitude_base_de_lancamento} rad, " \
               f"Longitude inicial: {self.longitude_base_de_lancamento} rad, " \
               f"Comprimento do trilho: {self.comprimento_do_trilho} m"


class ConstrutorBaseDeLancamento:

    def __init__(self):
        self.base_de_lancamento = BaseDeLancamento()

    def com_altitude_base(self, h0):
        self.base_de_lancamento.altitude_base_de_lancamento = h0
        return self

    def com_latitude_inicial(self, delta0):
        self.base_de_lancamento.latitude_base_de_lancamento = -2.3267844 * np.pi / 180 if delta0 is None else delta0
        return self

    def com_longitude_inicial(self, lon0):
        self.base_de_lancamento.longitude_base_de_lancamento = -44.4111042 * np.pi / 180 if lon0 is None else lon0
        return self

    def com_comprimento_trilho(self, l_trilho):
        self.base_de_lancamento.comprimento_do_trilho = l_trilho
        return self

    def construir(self):
        return self.base_de_lancamento


# Exemplo de uso
builder = ConstrutorBaseDeLancamento()
parametros = (builder.com_altitude_base(0)
              .com_latitude_inicial(None)
              .com_longitude_inicial(None)
              .com_comprimento_trilho(10)
              .construir())

print(parametros)
