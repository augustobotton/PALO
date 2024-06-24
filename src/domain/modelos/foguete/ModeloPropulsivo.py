import numpy as np

from src.domain.modelos.foguete.modelo_propulsores import *

TEMPO_DE_IGNICAO_3_ESTAGIO = 1e10


class ModeloPropulsivo:
    def __init__(self, builder):
        """
        Inicializa a classe ModeloPropulsivo com base no builder fornecido.

        :param builder: Objeto do tipo ModeloPropulsivoBuilder que contém os parâmetros de construção.
        """
        self.impulso_especifico = builder.impulso_especifico
        self.tempos_de_separacao = None  #ts
        self.tempos_de_fim_de_queima = None  #tq
        self.tempos_de_ignicao = None  #ti
        self.massa_inicial_do_foguete = None
        self.distancia_radial_inicial = None
        self.duracao_da_segunda_queima_3_estagio = None
        self.velocidade_de_exaustao = None
        self.massa_propelente_estagios = np.array(builder.massa_propelente_estagios)
        self.massa_total_propelente_terceiro_estagio = builder.massa_total_propelente_terceiro_estagio
        self.duracao_de_queima_primeiro_estagio = builder.duracao_queima_estagios[0]
        self.duracao_de_queima_segundo_estagio = builder.duracao_queima_estagios[1]
        self.duracao_total_de_queima_do_terceiro_estagio = builder.duracao_queima_estagios[2]
        self.tempo_de_espera_separacao_1_2 = builder.tempo_espera_separacao[0]
        self.tempo_de_espera_separacao_2_3 = builder.tempo_espera_separacao[1]
        self.tempo_de_espera_separacao_3_4 = builder.tempo_espera_separacao[2]
        self.tempo_espera_ignicao_2_estagio = builder.tempo_espera_ignicao[0]
        self.tempo_espera_ignicao_3_estagio = builder.tempo_espera_ignicao[1]
        self.duracao_da_primeira_queima_3_estagio = builder.tempo_primeira_queima_terceiro_estagio
        self.massa_estrutural_por_estagio = np.array(builder.massa_estrutural_por_estagio)
        self.massa_de_carga_util = builder.massa_de_carga_util
        self.h0 = builder.h0
        self.planeta = builder.planeta
        self._calcular_tempos()

    def _calcular_tempos(self) -> None:
        """
        Calcula os tempos de ignição, queima e separação para os estágios do foguete.
        """
        self.duracao_da_segunda_queima_3_estagio = self.duracao_total_de_queima_do_terceiro_estagio - self.duracao_da_primeira_queima_3_estagio

        self.tempos_de_fim_de_queima = np.zeros(5)
        self.tempos_de_separacao = np.zeros(3)
        self.tempos_de_ignicao = np.zeros(5)

        self.tempos_de_ignicao[0] = 0
        self.tempos_de_fim_de_queima[0] = self.tempos_de_ignicao[0] + self.duracao_de_queima_primeiro_estagio
        self.tempos_de_separacao[0] = self.tempos_de_fim_de_queima[0] + self.tempo_de_espera_separacao_1_2
        self.tempos_de_ignicao[1] = self.tempos_de_separacao[0] + self.tempo_espera_ignicao_2_estagio
        self.tempos_de_fim_de_queima[1] = self.tempos_de_ignicao[1] + self.duracao_de_queima_segundo_estagio
        self.tempos_de_separacao[1] = self.tempos_de_fim_de_queima[1] + self.tempo_de_espera_separacao_2_3
        self.tempos_de_ignicao[2] = self.tempos_de_separacao[
                                        1] + self.tempo_espera_ignicao_3_estagio  #cuidar se o tempo de espera ta correto
        self.tempos_de_fim_de_queima[2] = self.tempos_de_ignicao[2] + self.duracao_da_primeira_queima_3_estagio
        self.tempos_de_ignicao[3] = TEMPO_DE_IGNICAO_3_ESTAGIO
        self.tempos_de_fim_de_queima[3] = self.tempos_de_ignicao[
                                              3] + self.duracao_da_segunda_queima_3_estagio  #atualizar?
        self.tempos_de_separacao[2] = self.tempos_de_fim_de_queima[3] + self.tempo_de_espera_separacao_3_4

        self.massa_propelente_estagios[2] = (
                self.massa_total_propelente_terceiro_estagio * self.duracao_da_primeira_queima_3_estagio /
                self.duracao_total_de_queima_do_terceiro_estagio
        )
        self.massa_propelente_estagios[3] = (
                self.massa_total_propelente_terceiro_estagio * self.duracao_da_segunda_queima_3_estagio /
                self.duracao_total_de_queima_do_terceiro_estagio
        )

        self.massa_inicial_do_foguete = np.sum(self.massa_propelente_estagios) + np.sum(
            self.massa_estrutural_por_estagio) + self.massa_de_carga_util

        self.distancia_radial_inicial = self.planeta.raio_equatorial + self.h0

    def __str__(self):
        return f"impulso_especifico: {self.impulso_especifico}\n" \
               f"tempos_de_separacao: {self.tempos_de_separacao}\n" \
               f"tempos_de_fim_de_queima: {self.tempos_de_fim_de_queima}\n" \
               f"tempos_de_ignicao: {self.tempos_de_ignicao}\n" \
               f"massa_inicial_do_foguete: {self.massa_inicial_do_foguete}\n" \
               f"distancia_radial_inicial: {self.distancia_radial_inicial}\n" \
               f"duracao_da_segunda_queima_3_estagio: {self.duracao_da_segunda_queima_3_estagio}\n" \
               f"velocidade_de_exaustao: {self.velocidade_de_exaustao}\n" \
               f"massa_propelente_estagios: {self.massa_propelente_estagios}\n" \
               f"massa_total_propelente_terceiro_estagio: {self.massa_total_propelente_terceiro_estagio}\n" \
               f"duracao_de_queima_primeiro_estagio: {self.duracao_de_queima_primeiro_estagio}\n" \
               f"duracao_de_queima_segundo_estagio: {self.duracao_de_queima_segundo_estagio}\n" \
               f"duracao_total_de_queima_do_terceiro_estagio: {self.duracao_total_de_queima_do_terceiro_estagio}\n" \
               f"tempo_de_espera_separacao_1_2: {self.tempo_de_espera_separacao_1_2}\n" \
               f"tempo_de_espera_separacao_2_3: {self.tempo_de_espera_separacao_2_3}\n" \
               f"tempo_de_espera_separacao_3_4: {self.tempo_de_espera_separacao_3_4}\n" \
               f"tempo_espera_ignicao_2_estagio: {self.tempo_espera_ignicao_2_estagio}\n" \
               f"tempo_espera_ignicao_3_estagio: {self.tempo_espera_ignicao_3_estagio}\n" \
               f"duracao_da_primeira_queima_3_estagio: {self.duracao_da_primeira_queima_3_estagio}\n" \

    def propulsao_n_estagios(self, t: float, X: np.array) -> np.array:
        """
        Função para cálculo dos parâmetros propulsivos em função do tempo
        Veículo de até 3 estágios com dupla ignição do terceiro estágio
        Entrada:
            t (s): tempo
            X: vetor de estado
            V = X[0] (m/s): módulo do vetor velocidade relativa com respeito ao planeta girante
            A = X[1] (rad): ângulo de azimute do vetor velocidade relativa com respeito ao eixo z (que aponta para o norte) do sistema uen.
            phi = X[2] (rad): ângulo de elevação do vetor velocidade relativa com respeito ao horizonte local (plano yz do referencial uen)
            r = X[3] (m): distância radial até o centro do planeta
            delta = X[4] (rad): latitude com respeito ao plano equatorial do planeta
            lon = X[5] (rad): longitude planetária
        Saídas:
            ft (N): força propulsiva
            m (kg): massa do foguete em função do tempo
            mu, epsl (rad): Ângulos de apontamento da tubeira
        Hipóteses:
            - Em cada estágio, é assumida taxa de queima contínua, ou seja, não há controle da queima e a mesma é assumida uniforme do início ao fim do propelente
            - A tração de cada estágio é assumida como um pulso retangular perfeito, ou seja, quando acionado, o propulsor vai de tração zero até a máxima, permanecendo nesse patamar constante. Ao fim da queima, a tração cai instantaneamente a zero
        """
        # Número de estágios
        N = len(self.tempos_de_ignicao)

        if N == 1:
            ft, m = propulsor_1_estagio(t, self.tempos_de_ignicao, self.tempos_de_fim_de_queima,
                                        self.tempos_de_separacao,
                                        self.impulso_especifico, self.massa_propelente_estagios,
                                        self.massa_estrutural_por_estagio,
                                        self.massa_inicial_do_foguete, self.planeta.gravidade)
        elif N == 2:
            ft, m = propulsor_2_estagios(t, self.tempos_de_ignicao, self.tempos_de_fim_de_queima,
                                         self.tempos_de_separacao,
                                         self.impulso_especifico, self.massa_propelente_estagios,
                                         self.massa_estrutural_por_estagio,
                                         self.massa_inicial_do_foguete, self.planeta.gravidade)
        elif N == 3:
            ft, m = propulsor_3_estagios(t, self.tempos_de_ignicao, self.tempos_de_fim_de_queima,
                                         self.tempos_de_separacao,
                                         self.impulso_especifico, self.massa_propelente_estagios,
                                         self.massa_estrutural_por_estagio,
                                         self.massa_inicial_do_foguete, self.planeta.gravidade)
        else:
            ft, m = propulsor_3_estagios_2ig(t, self.tempos_de_ignicao, self.tempos_de_fim_de_queima,
                                             self.tempos_de_separacao,
                                             self.impulso_especifico, self.massa_propelente_estagios,
                                             self.massa_estrutural_por_estagio, self.massa_inicial_do_foguete,
                                             self.massa_de_carga_util, self.planeta.gravidade)

        # Para altitudes acima de 200km, alinha o vetor de tração com a velocidade inercial ao invés da relativa
        # Desmembra o vetor de estado
        V, A, phi, r, delta = X[0], X[1], X[2], X[3], X[4]

        # Altitude
        h = r - self.planeta.raio_equatorial

        if h < 200e3:
            epsl = 0
            mu = 0  # Tração alinhada com a velocidade relativa
        else:
            # Vetor velocidade inercial
            _, phii, Ai = Vrel2Vine(V, phi, A, self.planeta.velocidade_inercial_de_rotacao, r, delta)

            # Ângulos propulsivos para que a tração seja alinhada com a velocidade inercial
            mu = np.arcsin(np.cos(A) * np.cos(phii) * np.sin(Ai) - np.sin(A) * np.cos(phii) * np.cos(Ai))
            epsl = -np.arctan2(
                -np.cos(phi) * np.sin(phii) + np.sin(phi) * np.sin(A) * np.cos(phii) * np.sin(Ai) + np.sin(
                    phi) * np.cos(A) * np.cos(phii) * np.cos(Ai),
                np.sin(phi) * np.sin(phii) + np.cos(phi) * np.sin(A) * np.cos(phii) * np.sin(Ai) + np.cos(
                    phi) * np.cos(A) * np.cos(phii) * np.cos(Ai))
        return ft, m, mu, epsl


class ConstrutorModeloPropulsivo:
    def __init__(self):
        """
        Inicializa a classe ModeloPropulsivoBuilder com valores padrão.
        """
        self.impulso_especifico = np.array([])
        self.massa_propelente_estagios = np.array([])
        self.massa_total_propelente_terceiro_estagio = 0.0
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

    def com_massa_total_propelente_terceiro_estagio(self,
                                                    massa_total_propelente_terceiro_estagio: float) -> 'ConstrutorModeloPropulsivo':
        """
        Define a massa de propelente do terceiro estágio.

        :param massa_total_propelente_terceiro_estagio: Massa de propelente do terceiro estágio.
        :return: Instância do builder.
        """
        self.massa_total_propelente_terceiro_estagio = massa_total_propelente_terceiro_estagio
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
