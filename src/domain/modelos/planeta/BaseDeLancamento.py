
class BaseDeLancamento:
    """
    Classe para criar e gerenciar uma base de lançamento com parâmetros configuráveis.

    """
    def __init__(self):
        """
        Inicializa uma nova instância da classe BaseDeLancamento com todos os parâmetros definidos como None.
        Deve ser usado 'ConstrutorBaseDeLancamento 'para configurar os parâmetros da base de lançamento.
        """
        self.altitude_base_de_lancamento = None  # Altitude da base de lançamento (em metros)
        self.latitude_base_de_lancamento = None  # Latitude inicial da base de lançamento (em radianos)
        self.longitude_base_de_lancamento = None  # Longitude inicial da base de lançamento (em radianos)
        self.comprimento_do_trilho = None  # Comprimento do trilho de lançamento (em metros)

    def __str__(self):
        """
        Retorna uma representação string da base de lançamento.

        Returns:
            str: Uma string que descreve a base de lançamento com suas coordenadas e comprimento do trilho.
        """
        return (f"Altitude da base de lançamento: {self.altitude_base_de_lancamento} m, "
                f"Latitude inicial: {self.latitude_base_de_lancamento} rad, "
                f"Longitude inicial: {self.longitude_base_de_lancamento} rad, "
                f"Comprimento do trilho: {self.comprimento_do_trilho} m")
