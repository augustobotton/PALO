import numpy as np


class ModeloFoguete:
    def __init__(self):
        # Parâmetros de massa estrutural e de carga útil
        self.massa_estrutural_por_estagio = np.array([7385, 1367, 59.69])
        self.massa_de_carga_util = 13

        # Parâmetros aerodinâmicos e ambientais
        fator_correcao = 1.28  # Fator de correção do arrasto
        self.area_secao_transversal_1_estagio = 4.6 * 5 / 3
        self.area_secao_transversal_2_estagio = 1.5
        self.area_secao_transcersal_3_estagio = 1.5
        self.area_secao_transversal_carga_util = 1.5

        self.comprimento_total_do_foguete = 7.33 + 7.1 + 6.28
        self.comprimento_sem_1_estagio = 7.1 + 6.28
        self.comprimento_sem_2_estagio = 6.28
        self.comprimento_carga_util = 1

        fator_de_correcao_2_estagio = (
                                              self.comprimento_sem_1_estagio / self.comprimento_total_do_foguete) * 0.5 + 0.5
        fator_de_correcao_3_estagio = (
                                              self.comprimento_sem_2_estagio / self.comprimento_total_do_foguete) * 0.5 + 0.5
        fator_correcao_carga_util = (
                                            self.comprimento_carga_util / self.comprimento_total_do_foguete) * 0.5 + 0.5

        areas_de_referencia_para_calculo_do_arrasto = np.array(
            [self.area_secao_transversal_1_estagio,
             self.area_secao_transversal_2_estagio * fator_de_correcao_2_estagio,
             self.area_secao_transcersal_3_estagio * fator_de_correcao_3_estagio,
             self.area_secao_transversal_carga_util * fator_correcao_carga_util]).reshape(-1, 1)

        comprimento_caracteristico = 1.5  # Comprimento característico - diâmetro dos estágios 2 e superiores

    def estudo_delta_v(self, planeta):
        # Estudo simplificado pela equação de foguete
        self.mpx = np.array(
            [self.massa_propelente_estagios_1_2[0], self.massa_propelente_estagios_1_2[1],
             self.massa_propelente_terceiro_estagio])
        self.ms_mpx_sum = self.massa_estrutural_por_estagio + self.mpx
        self.sigma = self.massa_estrutural_por_estagio / self.ms_mpx_sum

        self.massa_total_decolagem = self.massa_inicial_do_foguete
        self.massa_total_na_ignicacao_2_estagio = self.massa_estrutural_por_estagio[1] + self.mpx[1] + \
                                                  self.massa_estrutural_por_estagio[
                                                      2] + \
                                                  self.mpx[2] + self.massa_de_carga_util
        self.m03 = self.massa_estrutural_por_estagio[2] + self.mpx[2] + self.massa_de_carga_util

        self.lamb0 = self.massa_total_na_ignicacao_2_estagio / self.massa_total_decolagem
        self.lamb1 = self.m03 / self.massa_total_na_ignicacao_2_estagio
        self.lamb2 = self.massa_de_carga_util / self.m03
        self.lamb = np.array([self.lamb0, self.lamb1, self.lamb2])

        self.lambL = np.prod(self.lamb)

        self.velocidade_de_exaustao = planeta.gravidade * self.impulso_especico_por_estagio

        self.Dv = -np.sum(self.velocidade_de_exaustao * np.log(self.sigma + (1 - self.sigma) * self.lamb))

    def mostra_dados(self):
        # Mostra dados na tela
        print('Area de referencia do foguete com primeiro estagio (m^2):', self.Sr[0])
        print('Area de referencia do foguete com segundo estagio (m^2):', self.Sr[1])
        print('Area de referencia do foguete com terceiro estagio (m^2):', self.Sr[2])
        print('Area de referencia da carga util (m^2):', self.Sr[3])
        print('Massa inicial antes da queima do primeiro estagio - kg:', self.massa_total_decolagem)
        print('Massa inicial antes da queima do segundo estagio - kg:', self.massa_total_na_ignicacao_2_estagio)
        print('Massa inicial antes da queima do terceiro estagio - kg:', self.m03)
        print('Massa da carga util - kg:', self.massa_de_carga_util)
        print('Razoes estruturais:', self.sigma)
        print('Razoes de carga util:', self.lamb)
        print('Velocidades de exaustao - m/s:', self.velocidade_de_exaustao)
        print('Razao de carga útil total:', self.lambL)
        print('Impulso de velocidade total ideal - m/s:', self.Dv)
