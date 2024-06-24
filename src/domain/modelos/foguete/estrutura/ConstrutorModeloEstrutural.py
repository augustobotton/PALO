from src.domain.modelos.foguete.estrutura.ModeloEstrutural import ModeloEstrutural


class ConstrutorModeloEstrutural:
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

    def com_massa_estrutural_por_estagio(self, massa_estrutural_por_estagio: list):
        self.massa_estrutural_por_estagio = massa_estrutural_por_estagio
        return self

    def com_massa_de_carga_util(self, massa_de_carga_util: float):
        self.massa_de_carga_util = massa_de_carga_util
        return self

    def com_area_secao_transversal_1_estagio(self, area_secao_transversal_1_estagio: float):
        self.area_secao_transversal_1_estagio = area_secao_transversal_1_estagio
        return self

    def com_area_secao_transversal_2_estagio(self, area_secao_transversal_2_estagio: float):
        self.area_secao_transversal_2_estagio = area_secao_transversal_2_estagio
        return self

    def com_area_secao_transversal_3_estagio(self, area_secao_transversal_3_estagio: float):
        self.area_secao_transversal_3_estagio = area_secao_transversal_3_estagio
        return self

    def com_area_secao_transversal_carga_util(self, area_secao_transversal_carga_util: float):
        self.area_secao_transversal_carga_util = area_secao_transversal_carga_util
        return self

    def com_comprimento_total_do_foguete(self, comprimento_total_do_foguete: float):
        self.comprimento_total_do_foguete = comprimento_total_do_foguete
        return self

    def com_comprimento_sem_1_estagio(self, comprimento_sem_1_estagio: float):
        self.comprimento_sem_1_estagio = comprimento_sem_1_estagio
        return self

    def com_comprimento_sem_2_estagio(self, comprimento_sem_2_estagio: float):
        self.comprimento_sem_2_estagio = comprimento_sem_2_estagio
        return self

    def com_comprimento_carga_util(self, comprimento_carga_util: float):
        self.comprimento_carga_util = comprimento_carga_util
        return self

    def construir(self) -> ModeloEstrutural:
        """
        Constrói e retorna uma instância de ModeloEstrutural com os parâmetros configurados.

        :return: Instância de ModeloEstrutural.
        """
        return ModeloEstrutural(self.massa_estrutural_por_estagio, self.massa_de_carga_util,
                                self.area_secao_transversal_1_estagio, self.area_secao_transversal_2_estagio,
                                self.area_secao_transversal_3_estagio, self.area_secao_transversal_carga_util,
                                self.comprimento_total_do_foguete, self.comprimento_sem_1_estagio,
                                self.comprimento_sem_2_estagio, self.comprimento_carga_util)
