import numpy as np

from src.domain.modelos.foguete.aerodinamica.ModeloAerodinamico import ModeloAerodinamico
from src.domain.modelos.foguete.estrutura.ModeloEstrutural import ModeloEstrutural
from src.domain.modelos.foguete.propulsao.ModeloPropulsivo import ModeloPropulsivo


class ModeloFoguete:
    def __init__(self, modelo_propulsivo: ModeloPropulsivo, modelo_estrutural: ModeloEstrutural,
                 modelo_aerodinamico: ModeloAerodinamico):
        """
        Inicializa a classe ModeloFoguete com os modelos fornecidos.

        :param modelo_propulsivo: Instância de ModeloPropulsivo
        :param modelo_estrutural: Instância de ModeloEstrutural
        :param modelo_aerodinamico: Instância de ModeloAerodinamico
        """

        self.modelo_propulsivo = modelo_propulsivo
        self.modelo_estrutural = modelo_estrutural
        self.modelo_aerodinamico = modelo_aerodinamico
        self.massa_de_carga_util = modelo_estrutural.massa_de_carga_util
        self.massa_inicial_do_foguete = self.massa_inicial()
        self.modelo_propulsivo.massa_inicial_do_foguete = self.massa_inicial()
        self.impulso_especifico_por_estagio = modelo_propulsivo.impulso_especifico
        self.areas_referencia = self.modelo_estrutural.areas_de_referencia_para_calculo_do_arrasto

    def _estudo_delta_v(self, planeta):
        """
        Realiza um estudo simplificado pela equação de foguete.

        :param planeta: Instância do planeta com o atributo gravidade
        """
        massa_prop_por_estagio = np.array(
            [self.modelo_propulsivo.massa_propelente_estagios[0], self.modelo_propulsivo.massa_propelente_estagios[1],
             (self.modelo_propulsivo.massa_propelente_estagios[2] + self.modelo_propulsivo.massa_propelente_estagios[3])])
        razoes_estruturais = self.modelo_estrutural.massa_estrutural_por_estagio / (
                self.modelo_estrutural.massa_estrutural_por_estagio + massa_prop_por_estagio)

        massa_total_na_ignicacao_2_estagio = self.modelo_estrutural.massa_estrutural_por_estagio[1] + \
                                             self.modelo_propulsivo.massa_propelente_estagios[
                                                 1] + \
                                             self.modelo_estrutural.massa_estrutural_por_estagio[2] + \
                                             (self.modelo_propulsivo.massa_propelente_estagios[
                                                 2] +self.modelo_propulsivo.massa_propelente_estagios[3]) + self.massa_de_carga_util
        massa_total_na_ignicao_3_estagio = self.modelo_estrutural.massa_estrutural_por_estagio[2] + \
                                           (self.modelo_propulsivo.massa_propelente_estagios[
                                               2] +self.modelo_propulsivo.massa_propelente_estagios[3])+ self.massa_de_carga_util

        razao_carga_util_primeiro_estagio = massa_total_na_ignicacao_2_estagio / self.massa_inicial_do_foguete
        razao_carga_util_segundo_estagio = massa_total_na_ignicao_3_estagio / massa_total_na_ignicacao_2_estagio
        razao_carga_util_terceiro_estagio = self.massa_de_carga_util / massa_total_na_ignicao_3_estagio
        razoes_carga_util = np.array(
            [razao_carga_util_primeiro_estagio, razao_carga_util_segundo_estagio, razao_carga_util_terceiro_estagio])

        razao_carga_util_total = np.prod(razoes_carga_util)

        velocidade_de_exaustao = planeta.gravidade * self.impulso_especifico_por_estagio

        delta_v_ideal = -np.sum(
            velocidade_de_exaustao * np.log(razoes_estruturais + (1 - razoes_estruturais) * razoes_carga_util))

        return delta_v_ideal, razoes_carga_util, razoes_estruturais, velocidade_de_exaustao, razao_carga_util_total,\
                                        massa_total_na_ignicacao_2_estagio, massa_total_na_ignicao_3_estagio

    def massa_inicial(self):
        return np.sum(self.modelo_propulsivo.massa_propelente_estagios) + np.sum(
            self.modelo_estrutural.massa_estrutural_por_estagio) + self.massa_de_carga_util

    def mostra_dados(self):
        """
        Exibe os dados do foguete na tela.
        """

        delta_v_ideal, razoes_carga_util, razoes_estruturais, velocidade_de_exaustao, razao_carga_util_total, \
                  massa_total_na_ignicacao_2_estagio, massa_total_na_ignicao_3_estagio = self._estudo_delta_v(
            self.modelo_propulsivo.planeta)
        print(" ")

        print('Estudo de delta-v ideal do foguete:')
        print('Área de referência do foguete com primeiro estágio (m^2):', self.areas_referencia[0])
        print('Área de referência do foguete com segundo estágio (m^2):', self.areas_referencia[1])
        print('Área de referência do foguete com terceiro estágio (m^2):', self.areas_referencia[2])
        print('Área de referência da carga útil (m^2):', self.areas_referencia[3])
        print('Massa inicial antes da queima do primeiro estágio (kg):', self.massa_inicial_do_foguete)
        print('Massa inicial antes da queima do segundo estágio (kg):', massa_total_na_ignicacao_2_estagio)
        print('Massa inicial antes da queima do terceiro estágio (kg):', massa_total_na_ignicao_3_estagio)
        print('Massa da carga útil (kg):', self.massa_de_carga_util)
        print('Razões estruturais:', razoes_estruturais)
        print('Razões de carga útil:', razoes_carga_util)
        print('Velocidades de exaustão (m/s):', velocidade_de_exaustao)
        print('Razão de carga útil total:', razao_carga_util_total)
        print('Impulso de velocidade total ideal (m/s):', delta_v_ideal)
        print('')
