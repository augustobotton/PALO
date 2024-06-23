import numpy as np
from src.domain.utilidades_mecanica_orbital.orbitalUtils.Converte import Vrel2Vine

TEMPO_DE_IGNICAO_3_ESTAGIO = 1e10


class ModeloPropulsivo:
    def __init__(self, builder):
        """
        Inicializa a classe ModeloPropulsivo com base no builder fornecido.

        :param builder: Objeto do tipo ModeloPropulsivoBuilder que contém os parâmetros de construção.
        """
        self.impulso_especifico = None
        self.tempos_de_separacao = None
        self.tempos_de_fim_de_queima = None
        self.tempos_de_ignicao = None
        self.massa_inicial_do_foguete = None
        self.distancia_radial_inicial = None
        self.tempo_da_segunda_queima_3_estagio = None
        self.velocidade_de_exaustao = None
        self.impulso_especifico_por_estagio = np.array(builder.impulso_especifico)
        self.massa_propelente_estagios_1_2 = np.array(builder.massa_propelente_estagios)
        self.massa_propelente_terceiro_estagio = builder.massa_propelente_terceiro_estagio
        self.duracao_de_queima_primeiro_estagio = builder.duracao_queima_estagios[0]
        self.duracao_de_queima_segundo_estagio = builder.duracao_queima_estagios[1]
        self.duracao_total_de_queima_do_terceiro_estagio = builder.duracao_queima_estagios[2]
        self.tempo_de_espera_separacao_1_2 = builder.tempo_espera_separacao[0]
        self.tempo_de_espera_separacao_2_3 = builder.tempo_espera_separacao[1]
        self.tempo_de_espera_separacao_3_4 = builder.tempo_espera_separacao[2]
        self.tempo_espera_ignicao_2_estagio = builder.tempo_espera_ignicao[0]
        self.tempo_espera_ignicao_3_estagio = builder.tempo_espera_ignicao[1]
        self.tempo_da_primeira_queima_3_estagio = builder.tempo_primeira_queima_terceiro_estagio
        self.massa_estrutural_por_estagio = np.array(builder.massa_estrutural_por_estagio)
        self.massa_de_carga_util = builder.massa_de_carga_util
        self.h0 = builder.h0
        self.planeta = builder.planeta
        self._calcular_tempos()

    def _calcular_tempos(self) -> None:
        """
        Calcula os tempos de ignição, queima e separação para os estágios do foguete.
        """
        self.tempo_da_segunda_queima_3_estagio = self.duracao_total_de_queima_do_terceiro_estagio - self.tempo_da_primeira_queima_3_estagio

        self.tempos_de_fim_de_queima = np.zeros(4)
        self.tempos_de_separacao = np.zeros(4)
        self.tempos_de_ignicao = np.zeros(4)

        self.tempos_de_fim_de_queima[0] = self.tempos_de_ignicao[0] + self.duracao_de_queima_primeiro_estagio
        self.tempos_de_separacao[0] = self.tempos_de_fim_de_queima[0] + self.tempo_de_espera_separacao_1_2
        self.tempos_de_ignicao[1] = self.tempos_de_separacao[0] + self.tempo_espera_ignicao_2_estagio
        self.tempos_de_fim_de_queima[1] = self.tempos_de_ignicao[1] + self.duracao_de_queima_segundo_estagio
        self.tempos_de_separacao[1] = self.tempos_de_fim_de_queima[1] + self.tempo_de_espera_separacao_2_3
        self.tempos_de_ignicao[2] = self.tempos_de_separacao[1] + self.tempo_espera_ignicao_3_estagio
        self.tempos_de_fim_de_queima[2] = self.tempos_de_ignicao[2] + self.tempo_da_primeira_queima_3_estagio
        self.tempos_de_ignicao[3] = TEMPO_DE_IGNICAO_3_ESTAGIO
        self.tempos_de_fim_de_queima[3] = self.tempos_de_ignicao[3] + self.tempo_da_segunda_queima_3_estagio
        self.tempos_de_separacao[2] = self.tempos_de_fim_de_queima[3] + self.tempo_de_espera_separacao_3_4

        self.massa_propelente_estagios_1_2[2] = (
                self.massa_propelente_terceiro_estagio * self.tempo_da_primeira_queima_3_estagio /
                self.duracao_total_de_queima_do_terceiro_estagio
        )
        self.massa_propelente_estagios_1_2[3] = (
                self.massa_propelente_terceiro_estagio * self.tempo_da_segunda_queima_3_estagio /
                self.duracao_total_de_queima_do_terceiro_estagio
        )

        self.massa_inicial_do_foguete = np.sum(self.massa_propelente_estagios_1_2) + np.sum(
            self.massa_estrutural_por_estagio) + self.massa_de_carga_util
        # self.distancia_radial_inicial = terra.raio_equatorial_terrestre + self.h0

    def calcular_empuxo_massa(self, tempo: float, tempos_de_ignicao: np.ndarray, tempos_de_fim_de_queima: np.ndarray,
                              tempos_de_separacao: np.ndarray, impulso_especifico: np.ndarray,
                              massa_propelente: np.ndarray,
                              massa_estrutural: np.ndarray, massa_inicial: float, gravidade: float, estagio: int,
                              massa_carga_util: float = None) -> tuple:
        """
        Calcula o empuxo e a massa para um determinado estágio do foguete.

        :param tempo: Tempo atual.
        :param tempos_de_ignicao: Array com os tempos de ignição dos estágios.
        :param tempos_de_fim_de_queima: Array com os tempos de fim de queima dos estágios.
        :param tempos_de_separacao: Array com os tempos de separação dos estágios.
        :param impulso_especifico: Array com os impulsos específicos dos estágios.

        :param massa_propelente: Array com as massas de propelente dos estágios.
        :param massa_estrutural: Array com as massas estruturais dos estágios.
        :param massa_inicial: Massa inicial do foguete.
        :param gravidade: Valor da gravidade.
        :param estagio: Número do estágio atual.
        :param massa_carga_util: Massa da carga útil (opcional).
        :return: Tupla de empuxo e massa.
        """
       # estagio -= 1  # Ajuste para arrays baseados em 0
        if tempo <= tempos_de_ignicao[estagio]:
            massa = massa_inicial if estagio == 0 else self.calcular_massa_apos_separacao(
                massa_propelente, massa_estrutural,
                massa_inicial, estagio)
            empuxo = 0
        elif tempo <= tempos_de_fim_de_queima[estagio]:
            taxa_de_fluxo_de_massa = -massa_propelente[estagio] / (
                    tempos_de_fim_de_queima[estagio] - tempos_de_ignicao[estagio])
            massa = massa_inicial + taxa_de_fluxo_de_massa * (
                    tempo - tempos_de_ignicao[estagio]) if estagio == 0 else self.calcular_massa_apos_separacao(
                massa_propelente,
                massa_estrutural,
                massa_inicial, estagio) + taxa_de_fluxo_de_massa * (tempo - tempos_de_ignicao[estagio])
            empuxo = -gravidade * impulso_especifico[estagio] * taxa_de_fluxo_de_massa
        elif tempo <= tempos_de_separacao[estagio]:
            massa = self.calcular_massa_apos_separacao(
                massa_propelente, massa_estrutural,
                massa_inicial, estagio)
            empuxo = 0
        else:
            massa = self.calcular_massa_apos_separacao(massa_propelente, massa_estrutural, massa_inicial, estagio + 1)
            massa = massa_carga_util if massa_carga_util and estagio == len(tempos_de_ignicao) - 1 else massa
            empuxo = 0
        return empuxo, massa

    @staticmethod
    def calcular_massa_apos_separacao(massa_propelente: np.ndarray,
                                      massa_estrutural: np.ndarray, massa_inicial: float, estagio: int) -> float:
        """
        Calcula a massa do foguete após a separação de um determinado estágio.

        :param massa_propelente: Array com as massas de propelente dos estágios.
        :param massa_estrutural: Array com as massas estruturais dos estágios.
        :param massa_inicial: Massa inicial do foguete.
        :param estagio: Número do estágio (baseado em 1).
        :return: Massa após a separação do estágio.
        """
        massa = massa_inicial
        for i in range(estagio):
            massa -= massa_propelente[i] + massa_estrutural[i]
        return massa

    def propulsao_n_estagios(self, tempo: float, vetor_de_estados: tuple) -> tuple:
        """
        Calcula a propulsão para N estágios do foguete.

        :param tempo: Tempo atual.
        :param vetor_de_estados: Estado atual do foguete (velocidade, ângulo, etc.).
        :return: Tupla de empuxo, massa, mu, epsl.
        """
        n_estagios = len(self.tempos_de_ignicao)

        # Determina o estágio atual com base no tempo
        estagio_atual = sum([tempo > self.tempos_de_separacao[i] for i in range(n_estagios)])

        # Calcula empuxo e massa para o estágio atual
        empuxo, massa = self.calcular_empuxo_massa(tempo, self.tempos_de_ignicao, self.tempos_de_fim_de_queima,
                                                   self.tempos_de_separacao,
                                                   self.impulso_especifico, self.massa_propelente_estagios_1_2,
                                                   self.massa_estrutural_por_estagio,
                                                   self.massa_inicial_do_foguete, self.planeta.gravidade, estagio_atual,
                                                   self.massa_de_carga_util)

        # Cálculos controle
        velocidade, angulo, phi, raio, delta, long = vetor_de_estados
        altura = raio - self.planeta.raio_equatorial
        if altura < 200e3:
            epsl = 0
            mu = 0
        else:
            _, phii, angulo_i = Vrel2Vine(velocidade, phi, angulo, self.planeta.velocidade_inercial_de_rotacao, raio, delta)
            mu = np.arcsin(
                np.cos(angulo) * np.cos(phii) * np.sin(angulo_i) - np.sin(angulo) * np.cos(phii) * np.cos(angulo_i))
            epsl = -np.arctan2(
                -np.cos(phi) * np.sin(phii) + np.sin(phi) * np.sin(angulo) * np.cos(phii) * np.sin(angulo_i) +
                np.sin(phi) * np.cos(angulo) * np.cos(phii) * np.cos(angulo_i), np.sin(phi) * np.sin(phii) +
                np.cos(phi) * np.sin(angulo) * np.cos(phii) * np.sin(angulo_i) + np.cos(phi) * np.cos(angulo) * np.cos(
                    phii) * np.cos(angulo_i))

        return empuxo, massa, mu, epsl


class ConstrutorModeloPropulsivo:
    def __init__(self):
        """
        Inicializa a classe ModeloPropulsivoBuilder com valores padrão.
        """
        self.impulso_especifico = np.array([])
        self.massa_propelente_estagios = np.array([])
        self.massa_propelente_terceiro_estagio = 0.0
        self.duracao_queima_estagios = np.array([])
        self.tempo_espera_separacao = np.array([])
        self.tempo_espera_ignicao = np.array([])
        self.tempo_primeira_queima_terceiro_estagio = 0.0
        self.massa_estrutural_por_estagio = np.array([])
        self.massa_de_carga_util = 0.0
        self.planeta = None
        self.h0 = 0.0

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

    def com_massa_propelente_terceiro_estagio(self,
                                              massa_propelente_terceiro_estagio: float) -> 'ConstrutorModeloPropulsivo':
        """
        Define a massa de propelente do terceiro estágio.

        :param massa_propelente_terceiro_estagio: Massa de propelente do terceiro estágio.
        :return: Instância do builder.
        """
        self.massa_propelente_terceiro_estagio = massa_propelente_terceiro_estagio
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

    def com_massa_estrutural_por_estagio(self, massa_estrutural_por_estagio: list) -> 'ConstrutorModeloPropulsivo':
        """
        Define a massa estrutural por estágio.

        :param massa_estrutural_por_estagio: Lista de massas estruturais.
        :return: Instância do builder.
        """
        self.massa_estrutural_por_estagio = np.array(massa_estrutural_por_estagio)
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
