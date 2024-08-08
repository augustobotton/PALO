import numpy as np

from src.domain.modelos.foguete.propulsao.ModeloPropulsivo import ModeloPropulsivo


class ConstrutorModeloPropulsivo:
    def __init__(self):
        """
        Inicializa a classe ModeloPropulsivoBuilder com valores padrão.
        """
        self.impulso_especifico = np.array([])
        self.massa_propelente_estagios = np.array([])
        self.duracao_queima_estagios = np.array([])
        self.tempo_espera_separacao = np.array([])
        self.tempo_espera_ignicao = np.array([])
        self.tempo_primeira_queima_terceiro_estagio = 0.0
        self.massa_de_carga_util = 0.0
        self.planeta = None
        self.h0 = 0.0
        self.numero_motores_por_estagio = np.array([])

    def com_impulso_especifico(self, impulso_especifico: list) -> 'ConstrutorModeloPropulsivo':
        """
        Define o impulso específico dos estágios.

        :param impulso_especifico: Lista de impulsos específicos.
        :return: Instância do builder.
        """
        self.impulso_especifico = np.array(impulso_especifico)
        return self

    def com_massa_propelente_estagios(self, massa_propelente_estagios: list) -> 'ConstrutorModeloPropulsivo':
        """
        Define a massa de propelente dos estágios.

        :param massa_propelente_estagios: Lista de massas de propelente.
        :return: Instância do builder.
        """
        self.massa_propelente_estagios = np.array(massa_propelente_estagios)
        return self


    def com_duracao_queima_estagios(self, duracao_queima_estagios: list) -> 'ConstrutorModeloPropulsivo':
        """
        Define a duração da queima dos estágios.

        :param duracao_queima_estagios: Lista de durações de queima.
        :return: Instância do builder.
        """
        self.duracao_queima_estagios = np.array(duracao_queima_estagios)
        return self

    def com_tempo_espera_separacao(self, tempo_espera_separacao: list) -> 'ConstrutorModeloPropulsivo':
        """
        Define o tempo de espera para a separação dos estágios.

        :param tempo_espera_separacao: Lista de tempos de espera para separação.
        :return: Instância do builder.
        """
        self.tempo_espera_separacao = np.array(tempo_espera_separacao)
        return self

    def com_tempo_espera_ignicao(self, tempo_espera_ignicao: list) -> 'ConstrutorModeloPropulsivo':
        """
        Define o tempo de espera para a ignição dos estágios.

        :param tempo_espera_ignicao: Lista de tempos de espera para ignição.
        :return: Instância do builder.
        """
        self.tempo_espera_ignicao = np.array(tempo_espera_ignicao)
        return self

    def com_tempo_primeira_queima_terceiro_estagio(self,
                                                   tempo_primeira_queima_terceiro_estagio: float) -> 'ConstrutorModeloPropulsivo':
        """
        Define o tempo da primeira queima do terceiro estágio.

        :param tempo_primeira_queima_terceiro_estagio: Tempo da primeira queima do terceiro estágio.
        :return: Instância do builder.
        """
        self.tempo_primeira_queima_terceiro_estagio = tempo_primeira_queima_terceiro_estagio
        return self


    def com_massa_de_carga_util(self, massa_de_carga_util: float) -> 'ConstrutorModeloPropulsivo':
        """
        Define a massa da carga útil.

        :param massa_de_carga_util: Massa da carga útil.
        :return: Instância do builder.
        """
        self.massa_de_carga_util = massa_de_carga_util
        return self

    def com_h0(self, h0: float) -> 'ConstrutorModeloPropulsivo':
        """
        Define a altura inicial do foguete.

        :param h0: Altura inicial.
        :return: Instância do builder.
        """
        self.h0 = h0
        return self

    def com_planeta(self, planeta: str) -> 'ConstrutorModeloPropulsivo':
        """
        Define o planeta do foguete.

        :param planeta: Nome do planeta.
        :return: Instância do builder.
        """
        self.planeta = planeta
        return self


    def construir(self) -> ModeloPropulsivo:
        """
        Constrói uma instância de ModeloPropulsivo com os parâmetros definidos.

        :return: Instância de ModeloPropulsivo.
        """
        return ModeloPropulsivo(self)
