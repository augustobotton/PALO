import numpy as np

class Orbita:
    def __init__(self, semi_eixo_maior: float = None, excentricidade: float = None, inclinacao: float = None,
                 raan: float = None, arg_periastro: float = None, anomalia_verdadeira: float = None,
                 tempo_de_periastro: float = 0, parametro: float = None):
        """
        Inicializa a classe Orbita com os parâmetros orbitais básicos.

        :param semi_eixo_maior: Semieixo maior da órbita (em metros).
        :param excentricidade: Excentricidade da órbita (adimensional).
        :param inclinacao: Inclinação da órbita.
        :param raan: Longitude do nó ascendente.
        :param arg_periastro: Argumento do periastro.
        :param anomalia_verdadeira: Anomalia verdadeira.
        """
        self.semi_eixo_maior = semi_eixo_maior
        self.excentricidade = excentricidade
        self.inclinacao = inclinacao
        self.raan = raan
        self.arg_periastro = arg_periastro
        self.anomalia_verdadeira = anomalia_verdadeira
        self.tempo_de_periastro = tempo_de_periastro
        # TODO criar metodo para overide mu
        self.parametro = parametro
        self.mu = 398600.4418e9  # Constante gravitacional da Terra, em m^3/s^2. É o padrão

    def __repr__(self):
        return (f"Orbita(semi_eixo_maior={self.semi_eixo_maior}, excentricidade={self.excentricidade}, "
                f"inclinacao={self.inclinacao}, raan={self.raan}, arg_periastro={self.arg_periastro}, "
                f"anomalia_verdadeira={self.anomalia_verdadeira}, tempo_periastro={self.tempo_de_periastro})")

    def calcular_parametro_orbital(self) -> float:
        """
        Calcula o parâmetro orbital.

        :return: Parâmetro orbital.
        """
        return self.semi_eixo_maior * (1 - self.excentricidade ** 2)

    def _calcular_distancia_orbital(self, p: float) -> float:
        """
        Calcula a distância orbital.

        :param p: Parâmetro orbital.
        :return: Distância orbital.
        """
        return p / (1 + self.excentricidade * np.cos(self.anomalia_verdadeira))

    def _calcular_posicao_orbital(self, r: float) -> np.ndarray:
        """
        Calcula a posição orbital.

        :param r: Distância orbital.
        :return: Vetor de posição no plano orbital.
        """
        x_orbital = r * np.cos(self.anomalia_verdadeira)
        y_orbital = r * np.sin(self.anomalia_verdadeira)
        z_orbital = 0
        return np.array([x_orbital, y_orbital, z_orbital])

    def _calcular_velocidade_orbital(self, p: float) -> np.ndarray:
        """
        Calcula a velocidade orbital.

        :param p: Parâmetro orbital.
        :return: Vetor de velocidade no plano orbital.
        """
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
        p = self.calcular_parametro_orbital()
        r = self._calcular_distancia_orbital(p)

        posicao_orbital = self._calcular_posicao_orbital(r)
        velocidade_orbital = self._calcular_velocidade_orbital(p)

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
