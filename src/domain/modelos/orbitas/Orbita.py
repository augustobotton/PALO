import numpy as np

from src.domain.modelos.orbitas.utilidades.vetor_de_estados_para_coe import vetor_de_estados_para_coe
from src.domain.modelos.orbitas.utilidades.funcoes_conversao import matriz_rotacao_orbital_inercial


class Orbita:
    def __init__(self, semi_eixo_maior: float = None, excentricidade: float = None, inclinacao: float = None,
                 raan: float = None, argumento_periastro: float = None, anomalia_verdadeira: float = None,
                 tempo_de_periastro: float = 0, parametro: float = None, quantidade_momento_angular = None):
        """
        Inicializa a classe Orbita com os parâmetros orbitais básicos.
        A unidade de comprimento padrão é km.

        """
        self.semi_eixo_maior = semi_eixo_maior
        self.excentricidade = excentricidade
        self.inclinacao = inclinacao
        self.raan = raan
        self.arg_periastro = argumento_periastro
        self.anomalia_verdadeira = anomalia_verdadeira
        self.tempo_de_periastro = tempo_de_periastro
        self.parametro = parametro
        self.quantidade_momento_angular = quantidade_momento_angular
        self.mu = 398600.4418  # Constante gravitacional da Terra, em km^3/s^2.

    def retorna_parametros(self):
        return self.semi_eixo_maior, self.excentricidade, self.inclinacao, self.raan, self.arg_periastro, self.anomalia_verdadeira, self.quantidade_momento_angular

    def __repr__(self):
        return (f"Orbita(semi_eixo_maior={self.semi_eixo_maior}, excentricidade={self.excentricidade}, "
                f"inclinacao={self.inclinacao}, raan={self.raan}, arg_periastro={self.arg_periastro}, "
                f"anomalia_verdadeira={self.anomalia_verdadeira})")

    def define_mu(self, mu: float):
        """
        Define a constante gravitacional para outro corpo celeste.

        :param mu: Nova constante gravitacional (em km^3/s^2).
        """
        self.mu = mu

    def calcular_parametro_orbital(self) -> float:

        return self.semi_eixo_maior * (1 - self.excentricidade ** 2)

    def _calcular_distancia_orbital(self) -> float:

        return self.calcular_parametro_orbital() / (1 + self.excentricidade * np.cos(self.anomalia_verdadeira))

    def calcular_vetor_posicao_orbital(self) -> np.ndarray:

        r = self._calcular_distancia_orbital()
        x_orbital = r * np.cos(self.anomalia_verdadeira)
        y_orbital = r * np.sin(self.anomalia_verdadeira)
        z_orbital = 0
        return np.array([x_orbital, y_orbital, z_orbital])

    def calcular_vetor_velocidade_orbital(self) -> np.ndarray:

        p = self.calcular_parametro_orbital()
        h = np.sqrt(self.mu * p)
        vx_orbital = -self.mu / h * np.sin(self.anomalia_verdadeira)
        vy_orbital = self.mu / h * (self.excentricidade + np.cos(self.anomalia_verdadeira))
        vz_orbital = 0
        return np.array([vx_orbital, vy_orbital, vz_orbital])

    def calcular_vetor_estado(self) -> tuple[np.ndarray, np.ndarray]:

        posicao_orbital = self.calcular_vetor_posicao_orbital()
        velocidade_orbital = self.calcular_vetor_velocidade_orbital()
        # Matriz de rotação considerando órbitas circulares e equatoriais
        if self.excentricidade == 0:
            # Órbita circular
            self.arg_periastro = 0
        if self.inclinacao == 0 or self.inclinacao == 180:
            # Órbita equatorial
            self.raan = 0

        matriz_rotacao = matriz_rotacao_orbital_inercial(self.raan, self.arg_periastro, self.inclinacao)
        posicao_inercial = np.dot(matriz_rotacao, posicao_orbital)
        velocidade_inercial = np.dot(matriz_rotacao, velocidade_orbital)

        return posicao_inercial, velocidade_inercial

    def calcula_apoastro(self):

        apoastro = self.semi_eixo_maior * (1 + self.excentricidade)
        return apoastro

    def calcula_periastro(self):

        periastro = self.semi_eixo_maior * (1 - self.excentricidade)
        return periastro

    @classmethod
    def circular(cls, semi_eixo_maior: float, inclinacao: float, ):
        """
        Cria uma instância da classe Orbita para uma órbita circular.

        """
        excentricidade = 0
        argumento_periastro = 0
        raan = 0
        anomalia_verdadeira = 0
        tempo_de_periastro = 0
        return cls(semi_eixo_maior=semi_eixo_maior, excentricidade=excentricidade, inclinacao=inclinacao,
                   raan=raan, argumento_periastro=argumento_periastro, anomalia_verdadeira=anomalia_verdadeira,
                   tempo_de_periastro=tempo_de_periastro)

    @classmethod
    def criar_pelo_vetor_de_estado(cls, posicao, velocidade, mu):
        eixo, exc, inclinacao, raan, argumento_de_periastro, anomalia_verdadeira, quantidade_momento_angular = vetor_de_estados_para_coe(posicao, velocidade, mu)
        return cls(semi_eixo_maior=eixo, excentricidade=exc, inclinacao=inclinacao, raan=raan,
                   argumento_periastro=argumento_de_periastro, anomalia_verdadeira=anomalia_verdadeira,
                   quantidade_momento_angular=quantidade_momento_angular)
