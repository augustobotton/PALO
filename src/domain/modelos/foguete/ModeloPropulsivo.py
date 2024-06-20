import numpy as np
from src.domain.modelos.planeta.ModeloPlaneta import terra
from src.domain.utilidades_mecanica_orbital.orbitalUtils.Converte import Vrel2Vine


class ModeloPropulsivo:
    def __init__(self, builder):
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

    def calcular_tempos(self) -> None:
        """Calcula os tempos de ignição, queima e separação para os estágios do foguete."""
        self.tempo_da_segunda_queima_3_estagio = self.duracao_total_de_queima_do_terceiro_estagio - self.tempo_da_primeira_queima_3_estagio

        tempos_de_ignicao = np.array([0, 0, 0, 0])
        tempos_de_fim_de_queima = np.array([0, 0, 0, 0])
        tempos_de_separacao = np.array([0, 0, 0, 0])

        tempos_de_fim_de_queima[0] = tempos_de_ignicao[0] + self.duracao_de_queima_primeiro_estagio
        tempos_de_separacao[0] = tempos_de_fim_de_queima[0] + self.tempo_de_espera_separacao_1_2
        tempos_de_ignicao[1] = tempos_de_separacao[0] + self.tempo_espera_ignicao_2_estagio
        tempos_de_fim_de_queima[1] = tempos_de_ignicao[1] + self.duracao_de_queima_segundo_estagio
        tempos_de_separacao[1] = tempos_de_fim_de_queima[1] + self.tempo_de_espera_separacao_2_3
        tempos_de_ignicao[2] = tempos_de_separacao[1] + self.tempo_espera_ignicao_3_estagio
        tempos_de_fim_de_queima[2] = tempos_de_ignicao[2] + self.tempo_da_primeira_queima_3_estagio
        tempos_de_ignicao[3] = 1e10
        tempos_de_fim_de_queima[3] = tempos_de_ignicao[3] + self.tempo_da_segunda_queima_3_estagio
        tempos_de_separacao[2] = tempos_de_fim_de_queima[3] + self.tempo_de_espera_separacao_3_4

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
        self.distancia_radial_inicial = terra.raio_equatorial_terrestre + self.h0

    def calcular_empuxo_massa(self, tempo: float, tempos_de_ignicao: np.ndarray, tempos_de_fim_de_queima: np.ndarray,
                              tempos_de_separacao: np.ndarray, impulso_especifico: np.ndarray, massa_propelente: np.ndarray,
                              massa_estrutural: np.ndarray, massa_inicial: float, gravidade: float, estagio: int,
                              massa_carga_util: float = None) -> tuple:
        """
        Calcula o empuxo e a massa para um determinado estágio do foguete.
        :param tempo: Tempo atual
        :return: Tupla de empuxo e massa
        """
        estagio -= 1  # Ajuste para arrays baseados em 0
        if tempo <= tempos_de_ignicao[estagio]:
            massa = massa_inicial if estagio == 0 else self.calcular_massa_apos_separacao(
                tempos_de_ignicao, tempos_de_fim_de_queima, tempos_de_separacao, massa_propelente, massa_estrutural, massa_inicial, estagio)
            empuxo = 0
        elif tempo <= tempos_de_fim_de_queima[estagio]:
            taxa_de_fluxo_de_massa = -massa_propelente[estagio] / (tempos_de_fim_de_queima[estagio] - tempos_de_ignicao[estagio])
            massa = massa_inicial + taxa_de_fluxo_de_massa * (tempo - tempos_de_ignicao[estagio]) if estagio == 0 else self.calcular_massa_apos_separacao(
                tempos_de_ignicao, tempos_de_fim_de_queima, tempos_de_separacao, massa_propelente, massa_estrutural, massa_inicial, estagio) + taxa_de_fluxo_de_massa * (tempo - tempos_de_ignicao[estagio])
            empuxo = -gravidade * impulso_especifico[estagio] * taxa_de_fluxo_de_massa
        elif tempo <= tempos_de_separacao[estagio]:
            massa = self.calcular_massa_apos_separacao(
                tempos_de_ignicao, tempos_de_fim_de_queima, tempos_de_separacao, massa_propelente, massa_estrutural, massa_inicial, estagio)
            empuxo = 0
        else:
            massa = self.calcular_massa_apos_separacao(
                tempos_de_ignicao, tempos_de_fim_de_queima, tempos_de_separacao, massa_propelente, massa_estrutural, massa_inicial, estagio + 1)
            massa = massa_carga_util if massa_carga_util and estagio == len(tempos_de_ignicao) - 1 else massa
            empuxo = 0
        return empuxo, massa

    def calcular_massa_apos_separacao(self, tempos_de_ignicao: np.ndarray, tempos_de_fim_de_queima: np.ndarray,
                                      tempos_de_separacao: np.ndarray, massa_propelente: np.ndarray,
                                      massa_estrutural: np.ndarray, massa_inicial: float, estagio: int) -> float:
        """
        Calcula a massa do foguete após a separação de um determinado estágio.
        :param estagio: Número do estágio (baseado em 1)
        :return: Massa após a separação do estágio
        """
        massa = massa_inicial
        for i in range(estagio):
            massa -= massa_propelente[i] + massa_estrutural[i]
        return massa

    def propulsao_n_estagios(self, tempo: float, vetor_de_estados: tuple, tempos_de_ignicao: np.ndarray,
                             tempos_de_fim_de_queima: np.ndarray, tempos_de_separacao: np.ndarray, impulso_especifico: np.ndarray,
                             massa_propelente: np.ndarray, massa_estrutural: np.ndarray, massa_inicial: float,
                             gravidade: float, velocidade_angular_terra: float, raio_terra: float,
                             massa_carga_util: float = None) -> tuple:
        """
        Calcula a propulsão para N estágios.

        :param tempo: Tempo atual
        :return: Tupla de empuxo, massa, mu, epsl
        """
        n_estagios = len(tempos_de_ignicao)

        # Determina o estágio atual com base no tempo
        estagio_atual = sum([tempo > tempos_de_separacao[i] for i in range(n_estagios)])

        # Calcula empuxo e massa para o estágio atual
        empuxo, massa = self.calcular_empuxo_massa(tempo, tempos_de_ignicao, tempos_de_fim_de_queima, tempos_de_separacao, impulso_especifico, massa_propelente, massa_estrutural, massa_inicial, gravidade, estagio_atual, massa_carga_util)

        # Cálculos adicionais
        velocidade, angulo, phi, raio, delta = vetor_de_estados
        altura = raio - raio_terra
        if altura < 200e3:
            epsl = 0
            mu = 0
        else:
            _, phii, angulo_i = Vrel2Vine(velocidade, phi, angulo, velocidade_angular_terra, raio, delta)
            mu = np.arcsin(np.cos(angulo) * np.cos(phii) * np.sin(angulo_i) - np.sin(angulo) * np.cos(phii) * np.cos(angulo_i))
            epsl = -np.arctan2(-np.cos(phi) * np.sin(phii) + np.sin(phi) * np.sin(angulo) * np.cos(phii) * np.sin(angulo_i) +
                               np.sin(phi) * np.cos(angulo) * np.cos(phii) * np.cos(angulo_i), np.sin(phi) * np.sin(phii) +
                               np.cos(phi) * np.sin(angulo) * np.cos(phii) * np.sin(angulo_i) + np.cos(phi) * np.cos(angulo) * np.cos(
                phii) * np.cos(angulo_i))

        return empuxo, massa, mu, epsl

class ModeloPropulsivoBuilder:
    def __init__(self):
        self.impulso_especifico = []
        self.massa_propelente_estagios = []
        self.massa_propelente_terceiro_estagio = 0.0
        self.duracao_queima_estagios = []
        self.tempo_espera_separacao = []
        self.tempo_espera_ignicao = []
        self.tempo_primeira_queima_terceiro_estagio = 0.0
        self.massa_estrutural_por_estagio = []
        self.massa_de_carga_util = 0.0
        self.h0 = 0.0

    def set_impulso_especifico(self, impulso_especifico: list):
        self.impulso_especifico = impulso_especifico
        return self

    def set_massa_propelente_estagios(self, massa_propelente_estagios: list):
        self.massa_propelente_estagios = massa_propelente_estagios
        return self

    def set_massa_propelente_terceiro_estagio(self, massa_propelente_terceiro_estagio: float):
        self.massa_propelente_terceiro_estagio = massa_propelente_terceiro_estagio
        return self

    def set_duracao_queima_estagios(self, duracao_queima_estagios: list):
        self.duracao_queima_estagios = duracao_queima_estagios
        return self

    def set_tempo_espera_separacao(self, tempo_espera_separacao: list):
        self.tempo_espera_separacao = tempo_espera_separacao
        return self

    def set_tempo_espera_ignicao(self, tempo_espera_ignicao: list):
        self.tempo_espera_ignicao = tempo_espera_ignicao
        return self

    def set_tempo_primeira_queima_terceiro_estagio(self, tempo_primeira_queima_terceiro_estagio: float):
        self.tempo_primeira_queima_terceiro_estagio = tempo_primeira_queima_terceiro_estagio
        return self

    def set_massa_estrutural_por_estagio(self, massa_estrutural_por_estagio: list):
        self.massa_estrutural_por_estagio = massa_estrutural_por_estagio
        return self

    def set_massa_de_carga_util(self, massa_de_carga_util: float):
        self.massa_de_carga_util = massa_de_carga_util
        return self

    def set_h0(self, h0: float):
        self.h0 = h0
        return self

    def build(self):
        return ModeloPropulsivo(self)
