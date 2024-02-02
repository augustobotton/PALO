import numpy as np


class ModeloEstrutural:
    def __init__(self, ):
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

        return areas_de_referencia_para_calculo_do_arrasto, comprimento_caracteristico, fator_correcao
