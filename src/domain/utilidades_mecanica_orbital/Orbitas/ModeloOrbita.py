import numpy as np

from src.domain.utilidades_mecanica_orbital.orbitalUtils.calculos_orbitais import determina_parametros_orbitais


class Orbita:
    def __init__(self, semi_eixo_maior: float = None, excentricidade: float = None, inclinacao: float = None,
                 raan: float = None, argumento_periastro: float = None, anomalia_verdadeira: float = None,
                 tempo_de_periastro: float = 0, parametro: float = None):
        """
        Inicializa a classe Orbita com os parâmetros orbitais básicos.

        :param semi_eixo_maior: Semieixo maior da órbita (em metros).
        :param excentricidade: Excentricidade da órbita (adimensional).
        :param inclinacao: Inclinação da órbita.
        :param raan: Longitude do nó ascendente.
        :param argumento_periastro: Argumento do periastro.
        :param anomalia_verdadeira: Anomalia verdadeira.
        """
        self.semi_eixo_maior = semi_eixo_maior
        self.excentricidade = excentricidade
        self.inclinacao = inclinacao
        self.raan = raan
        self.arg_periastro = argumento_periastro
        self.anomalia_verdadeira = anomalia_verdadeira
        self.tempo_de_periastro = tempo_de_periastro
        self.parametro = parametro
        self.mu = 398600.4418e9  # Constante gravitacional da Terra, em m^3/s^2. É o padrão

    def __repr__(self):
        return (f"Orbita(semi_eixo_maior={self.semi_eixo_maior}, excentricidade={self.excentricidade}, "
                f"inclinacao={self.inclinacao}, raan={self.raan}, arg_periastro={self.arg_periastro}, "
                f"anomalia_verdadeira={self.anomalia_verdadeira}, tempo_periastro={self.tempo_de_periastro})")

    def define_mu(self, mu: float):
        """
        Define a constante gravitacional para outro corpo celeste.

        :param mu: Nova constante gravitacional (em m^3/s^2).
        """
        self.mu = mu

    def calcular_parametro_orbital(self) -> float:
        """
        Calcula o parâmetro orbital.

        :return: Parâmetro orbital.
        """
        return self.semi_eixo_maior * (1 - self.excentricidade ** 2)

    def _calcular_distancia_orbital(self) -> float:
        """
        Calcula a distância orbital.

        :param p: Parâmetro orbital.
        :return: Distância orbital.
        """
        return self.calcular_parametro_orbital() / (1 + self.excentricidade * np.cos(self.anomalia_verdadeira))

    def calcular_vetor_posicao_orbital(self) -> np.ndarray:
        """
        Calcula a posição orbital.

        :param r: Distância orbital.
        :return: Vetor de posição no plano orbital.
        """
        r = self._calcular_distancia_orbital()
        x_orbital = r * np.cos(self.anomalia_verdadeira)
        y_orbital = r * np.sin(self.anomalia_verdadeira)
        z_orbital = 0
        return np.array([x_orbital, y_orbital, z_orbital])

    def calcular_vetor_velocidade_orbital(self) -> np.ndarray:
        """
        Calcula a velocidade orbital.

        :param p: Parâmetro orbital.
        :return: Vetor de velocidade no plano orbital.
        """
        p = self.calcular_parametro_orbital()
        h = np.sqrt(self.mu * p)
        vx_orbital = -self.mu / h * np.sin(self.anomalia_verdadeira)
        vy_orbital = self.mu / h * (self.excentricidade + np.cos(self.anomalia_verdadeira))
        vz_orbital = 0
        return np.array([vx_orbital, vy_orbital, vz_orbital])

    def calcular_vetor_estado(self) -> tuple[np.ndarray, np.ndarray]:
        """
        Calcula o vetor de estado (posição e velocidade) a partir dos elementos orbitais.

        :return: Vetores de posição e velocidade no sistema de coordenadas inerciais.
        """

        posicao_orbital = self.calcular_vetor_posicao_orbital()
        velocidade_orbital = self.calcular_vetor_velocidade_orbital()

        # Matriz de rotação considerando órbitas circulares e equatoriais
        if self.excentricidade == 0:
            # Órbita circular
            self.arg_periastro = 0
        if self.inclinacao == 0 or self.inclinacao == 180:
            # Órbita equatorial
            self.raan = 0

        matriz_rotacao = self._matriz_rotacao_orbital()

        posicao_inercial = np.dot(matriz_rotacao, posicao_orbital)
        velocidade_inercial = np.dot(matriz_rotacao, velocidade_orbital)

        return posicao_inercial, velocidade_inercial

    def _matriz_rotacao_orbital(self) -> np.ndarray:
        """
        Calcula a matriz de rotação do plano orbital para o sistema inercial.

        :return: Matriz de rotação.
        """
        cos_raan = np.cos(self.raan)
        sin_raan = np.sin(self.raan)
        cos_arg_periastro = np.cos(self.arg_periastro)
        sin_arg_periastro = np.sin(self.arg_periastro)
        cos_inclinacao = np.cos(self.inclinacao)
        sin_inclinacao = np.sin(self.inclinacao)

        matriz_rotacao = np.array([
            [cos_raan * cos_arg_periastro - sin_raan * sin_arg_periastro * cos_inclinacao,
             -cos_raan * sin_arg_periastro - sin_raan * cos_arg_periastro * cos_inclinacao,
             sin_raan * sin_inclinacao],
            [sin_raan * cos_arg_periastro + cos_raan * sin_arg_periastro * cos_inclinacao,
             -sin_raan * sin_arg_periastro + cos_raan * cos_arg_periastro * cos_inclinacao,
             -cos_raan * sin_inclinacao],
            [sin_arg_periastro * sin_inclinacao,
             cos_arg_periastro * sin_inclinacao,
             cos_inclinacao]
        ])

        return matriz_rotacao

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
        :param semi_eixo_maior: Semieixo maior da órbita (em metros).
        :param inclinacao: Inclinação da órbita.
        :param raan: Longitude do nó ascendente.
        :param anomalia_verdadeira: Anomalia verdadeira.
        :param tempo_de_periastro: Tempo do periastro (em segundos).
        :return: Instância da classe Orbita.
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
    def criar_pelo_vetor_de_estado(cls, tempo_observacao, mu, posicao, velocidade):
        eixo, exc, inclinacao, raan, argumento_de_periastro, anomalia_verdadeira, tempo_periastro = determina_parametros_orbitais(tempo_observacao, mu, posicao, velocidade)
        return cls (semi_eixo_maior=eixo, excentricidade=exc, inclinacao=inclinacao, raan=raan, argumento_periastro= argumento_de_periastro, anomalia_verdadeira=anomalia_verdadeira, tempo_de_periastro= tempo_periastro)