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






