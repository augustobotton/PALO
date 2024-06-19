import numpy as np
from src.domain.modelos.foguete.ModeloAerodinamico import ModeloAerodinamico
from src.domain.modelos.foguete.ModeloEstrutural import ModeloEstrutural
from src.domain.modelos.foguete.ModeloPropulsivo import ModeloPropulsivo

class ModeloFoguete:
    def __init__(self, modelo_propulsivo: ModeloPropulsivo, modelo_estrutural: ModeloEstrutural, modelo_aerodinamico: ModeloAerodinamico):
        """
        Inicializa a classe ModeloFoguete com os modelos fornecidos.

        :param modelo_propulsivo: Instância de ModeloPropulsivo
        :param modelo_estrutural: Instância de ModeloEstrutural
        :param modelo_aerodinamico: Instância de ModeloAerodinamico
        """
        self.modelo_propulsivo = modelo_propulsivo
        self.modelo_estrutural = modelo_estrutural
        self.modelo_aerodinamico = modelo_aerodinamico
        self.massa_propelente_estagios_1_2 = modelo_propulsivo.massa_propelente_estagios_1_2
        self.massa_propelente_terceiro_estagio = modelo_propulsivo.massa_propelente_terceiro_estagio
        self.massa_estrutural_por_estagio = modelo_estrutural.massa_estrutural_por_estagio
        self.massa_de_carga_util = modelo_estrutural.massa_de_carga_util
        self.massa_inicial_do_foguete = modelo_propulsivo.massa_inicial_do_foguete
        self.impulso_especifico_por_estagio = modelo_propulsivo.impulso_especifico_por_estagio
        self.Sr = [modelo_estrutural.area_secao_transversal_1_estagio,
                   modelo_estrutural.area_secao_transversal_2_estagio,
                   modelo_estrutural.area_secao_transversal_3_estagio,
                   modelo_estrutural.area_secao_transversal_carga_util]

    def estudo_delta_v(self, planeta):
        """
        Realiza um estudo simplificado pela equação de foguete.

        :param planeta: Instância do planeta com o atributo gravidade
        """
        mpx = np.array([self.massa_propelente_estagios_1_2[0], self.massa_propelente_estagios_1_2[1], self.massa_propelente_terceiro_estagio])
        ms_mpx_sum = self.massa_estrutural_por_estagio + mpx
        sigma = self.massa_estrutural_por_estagio / ms_mpx_sum

        massa_total_na_ignicacao_2_estagio = self.massa_estrutural_por_estagio[1] + mpx[1] + self.massa_estrutural_por_estagio[2] + mpx[2] + self.massa_de_carga_util
        m03 = self.massa_estrutural_por_estagio[2] + mpx[2] + self.massa_de_carga_util

        lamb0 = massa_total_na_ignicacao_2_estagio / self.massa_inicial_do_foguete
        lamb1 = m03 / massa_total_na_ignicacao_2_estagio
        lamb2 = self.massa_de_carga_util / m03
        lamb = np.array([lamb0, lamb1, lamb2])

        lambL = np.prod(lamb)

        velocidade_de_exaustao = planeta.gravidade * self.impulso_especifico_por_estagio

        Dv = -np.sum(velocidade_de_exaustao * np.log(sigma + (1 - sigma) * lamb))

        # Armazena os resultados como atributos da classe
        self.mpx = mpx
        self.ms_mpx_sum = ms_mpx_sum
        self.sigma = sigma
        self.massa_total_na_ignicacao_2_estagio = massa_total_na_ignicacao_2_estagio
        self.m03 = m03
        self.lamb0 = lamb0
        self.lamb1 = lamb1
        self.lamb2 = lamb2
        self.lamb = lamb
        self.lambL = lambL
        self.velocidade_de_exaustao = velocidade_de_exaustao
        self.Dv = Dv

    def mostra_dados(self):
        """
        Exibe os dados do foguete na tela.
        """
        print('Área de referência do foguete com primeiro estágio (m^2):', self.Sr[0])
        print('Área de referência do foguete com segundo estágio (m^2):', self.Sr[1])
        print('Área de referência do foguete com terceiro estágio (m^2):', self.Sr[2])
        print('Área de referência da carga útil (m^2):', self.Sr[3])
        print('Massa inicial antes da queima do primeiro estágio (kg):', self.massa_inicial_do_foguete)
        print('Massa inicial antes da queima do segundo estágio (kg):', self.massa_total_na_ignicacao_2_estagio)
        print('Massa inicial antes da queima do terceiro estágio (kg):', self.m03)
        print('Massa da carga útil (kg):', self.massa_de_carga_util)
        print('Razões estruturais:', self.sigma)
        print('Razões de carga útil:', self.lamb)
        print('Velocidades de exaustão (m/s):', self.velocidade_de_exaustao)
        print('Razão de carga útil total:', self.lambL)
        print('Impulso de velocidade total ideal (m/s):', self.Dv)


class ConstrutorDeFoguete:
    def __init__(self):
        self.modelo_propulsivo = None
        self.modelo_estrutural = None
        self.modelo_aerodinamico = None

    def com_modelo_propulsivo(self, modelo_propulsivo: ModeloPropulsivo):
        self.modelo_propulsivo = modelo_propulsivo
        return self

    def com_modelo_estrutural(self, modelo_estrutural: ModeloEstrutural):
        self.modelo_estrutural = modelo_estrutural
        return self

    def com_modelo_aerodinamico(self, modelo_aerodinamico: ModeloAerodinamico):
        self.modelo_aerodinamico = modelo_aerodinamico
        return self

    def construir(self):
        if not self.modelo_propulsivo or not self.modelo_estrutural or not self.modelo_aerodinamico:
            raise ValueError("Todos os modelos (propulsivo, estrutural, aerodinâmico) devem ser fornecidos.")
        return ModeloFoguete(self.modelo_propulsivo, self.modelo_estrutural, self.modelo_aerodinamico)
