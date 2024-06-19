import numpy as np


class BaseDeLancamento:
    def __init__(self):
        # Inicializa os atributos da base de lançamento
        self.altitude_base_de_lancamento = None  # Altitude da base de lançamento (em metros)
        self.latitude_base_de_lancamento = None  # Latitude inicial da base de lançamento (em radianos)
        self.longitude_base_de_lancamento = None  # Longitude inicial da base de lançamento (em radianos)
        self.comprimento_do_trilho = None  # Comprimento do trilho de lançamento (em metros)

    def __str__(self):
        # Retorna uma representação textual da base de lançamento
        return f"Altitude da base de lançamento: {self.altitude_base_de_lancamento} m, " \
               f"Latitude inicial: {self.latitude_base_de_lancamento} rad, " \
               f"Longitude inicial: {self.longitude_base_de_lancamento} rad, " \
               f"Comprimento do trilho: {self.comprimento_do_trilho} m"


class ConstrutorBaseDeLancamento:
    DEFAULT_LATITUDE = -2.3267844 * np.pi / 180  # Latitude padrão em radianos
    DEFAULT_LONGITUDE = -44.4111042 * np.pi / 180  # Longitude padrão em radianos

    def __init__(self):
        # Inicializa uma nova instância de BaseDeLancamento
        self.base_de_lancamento = BaseDeLancamento()

    def com_altitude_base(self, h0: float):
        # Configura a altitude da base de lançamento
        if h0 < 0:
            raise ValueError("A altitude da base de lançamento não pode ser negativa.")
        self.base_de_lancamento.altitude_base_de_lancamento = h0
        return self

    def com_latitude_inicial(self, delta0: float = None):
        # Configura a latitude inicial da base de lançamento
        self.base_de_lancamento.latitude_base_de_lancamento = self.DEFAULT_LATITUDE if delta0 is None else delta0
        return self

    def com_longitude_inicial(self, lon0: float = None):
        # Configura a longitude inicial da base de lançamento
        self.base_de_lancamento.longitude_base_de_lancamento = self.DEFAULT_LONGITUDE if lon0 is None else lon0
        return self

    def com_comprimento_trilho(self, l_trilho: float):
        # Configura o comprimento do trilho de lançamento
        if l_trilho <= 0:
            raise ValueError("O comprimento do trilho deve ser maior que zero.")
        self.base_de_lancamento.comprimento_do_trilho = l_trilho
        return self

    def construir(self) -> BaseDeLancamento:
        # Retorna a instância de BaseDeLancamento configurada
        return self.base_de_lancamento



