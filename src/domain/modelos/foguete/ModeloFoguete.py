import numpy as np

from src.domain.modelos.foguete.ModeloAerodinamico import ModeloAerodinamico
from src.domain.modelos.foguete.ModeloEstrutural import ModeloEstrutural
from src.domain.modelos.foguete.ModeloPropulsivo import ModeloPropulsivo


class ModeloFoguete:
    def __init__(self):
        self.modelo_propulsivo = ModeloPropulsivo()
        self.modelo_propulsivo = ModeloEstrutural()
        self.modelo_propulsivo = ModeloAerodinamico(self.modelo_propulsivo)

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
