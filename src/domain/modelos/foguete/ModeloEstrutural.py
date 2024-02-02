import numpy as np


class ModeloEstrutural:
    def __init__(self, massa_estrutural_por_estagio, massa_de_carga_util, area_secao_transversal_1_estagio,
                 area_secao_transversal_2_estagio, area_secao_transversal_3_estagio, area_secao_transversal_carga_util,
                 comprimento_total_do_foguete, comprimento_sem_1_estagio, comprimento_sem_2_estagio,
                 comprimento_carga_util):
        self.massa_estrutural_por_estagio = massa_estrutural_por_estagio
        self.massa_de_carga_util = massa_de_carga_util
        self.area_secao_transversal_1_estagio = area_secao_transversal_1_estagio
        self.area_secao_transversal_2_estagio = area_secao_transversal_2_estagio
        self.area_secao_transcersal_3_estagio = area_secao_transversal_3_estagio
        self.area_secao_transversal_carga_util = area_secao_transversal_carga_util
        self.comprimento_total_do_foguete = comprimento_total_do_foguete
        self.comprimento_sem_1_estagio = comprimento_sem_1_estagio
        self.comprimento_sem_2_estagio = comprimento_sem_2_estagio
        self.comprimento_carga_util = comprimento_carga_util

    def calcula(self):
        fator_correcao = 1.28  # Fator de correção do arrasto

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


class ModeloEstruturalBuilder:
    def __init__(self):
        self.massa_estrutural_por_estagio = None
        self.massa_de_carga_util = None
        self.area_secao_transversal_1_estagio = None
        self.area_secao_transversal_2_estagio = None
        self.area_secao_transversal_3_estagio = None
        self.area_secao_transversal_carga_util = None
        self.comprimento_total_do_foguete = None
        self.comprimento_sem_1_estagio = None
        self.comprimento_sem_2_estagio = None
        self.comprimento_carga_util = None

    def with_massa_estrutural_por_estagio(self, massa_estrutural_por_estagio):
        self.massa_estrutural_por_estagio = massa_estrutural_por_estagio
        return self

    def with_massa_de_carga_util(self, massa_de_carga_util):
        self.massa_de_carga_util = massa_de_carga_util
        return self

    def with_area_secao_transversal_1_estagio(self, area_secao_transversal_1_estagio):
        self.area_secao_transversal_1_estagio = area_secao_transversal_1_estagio
        return self

    def with_area_secao_transversal_2_estagio(self, area_secao_transversal_2_estagio):
        self.area_secao_transversal_2_estagio = area_secao_transversal_2_estagio
        return self

    def with_area_secao_transversal_3_estagio(self, area_secao_transversal_3_estagio):
        self.area_secao_transversal_3_estagio = area_secao_transversal_3_estagio
        return self

    def with_area_secao_transversal_carga_util(self, area_secao_transversal_carga_util):
        self.area_secao_transversal_carga_util = area_secao_transversal_carga_util
        return self

    def with_comprimento_total_do_foguete(self, comprimento_total_do_foguete):
        self.comprimento_total_do_foguete = comprimento_total_do_foguete
        return self

    def with_comprimento_sem_1_estagio(self, comprimento_sem_1_estagio):
        self.comprimento_sem_1_estagio = comprimento_sem_1_estagio
        return self

    def with_comprimento_sem_2_estagio(self, comprimento_sem_2_estagio):
        self.comprimento_sem_2_estagio = comprimento_sem_2_estagio
        return self

    def with_comprimento_carga_util(self, comprimento_carga_util):
        self.comprimento_carga_util = comprimento_carga_util
        return self

    def build(self):
        return ModeloEstrutural(self.massa_estrutural_por_estagio, self.massa_de_carga_util,
                                self.area_secao_transversal_1_estagio, self.area_secao_transversal_2_estagio,
                                self.area_secao_transversal_3_estagio, self.area_secao_transversal_carga_util,
                                self.comprimento_total_do_foguete, self.comprimento_sem_1_estagio,
                                self.comprimento_sem_2_estagio, self.comprimento_carga_util)
