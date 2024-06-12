import numpy as np
class Orbita:
    def __init__(self, semi_eixo_maior: float, excentricidade: float, inclinacao: float,
                 raan: float, arg_perigeu: float, anomalia_verdadeira: float):
        """
        Inicializa a classe Orbita com os parâmetros orbitais básicos.

        :param semi_eixo_maior: Semieixo maior da órbita (em metros).
        :param excentricidade: Excentricidade da órbita (adimensional).
        :param inclinacao: Inclinação da órbita (em graus).
        :param raan: Longitude do nó ascendente (em graus).
        :param arg_perigeu: Argumento do perigeu (em graus).
        :param anomalia_verdadeira: Anomalia verdadeira (em graus).
        """
        self.semi_eixo_maior = semi_eixo_maior
        self.excentricidade = excentricidade
        self.inclinacao = np.radians(inclinacao)
        self.raan = np.radians(raan)
        self.arg_perigeu = np.radians(arg_perigeu)
        self.anomalia_verdadeira = np.radians(anomalia_verdadeira)
        self.mu = 398600.4418e9  # Constante gravitacional da Terra, em m^3/s^2

    def calcular_vetor_estado(self) -> tuple[np.ndarray, np.ndarray]:
        """
        Calcula o vetor de estado (posição e velocidade) a partir dos elementos orbitais.

        :return: Vetores de posição e velocidade no sistema de coordenadas inerciais.
        """
        p = self._calcular_parametro_orbital()
        r = self._calcular_distancia_orbital(p)

        posicao_orbital = self._calcular_posicao_orbital(r)
        velocidade_orbital = self._calcular_velocidade_orbital(p)

        matriz_rotacao = UtilidadesOrbitais.matriz_rotacao_orbital_inercial(
            self.raan, self.arg_perigeu, self.inclinacao)

        posicao_inercial = np.dot(matriz_rotacao, posicao_orbital)
        velocidade_inercial = np.dot(matriz_rotacao, velocidade_orbital)

        return posicao_inercial, velocidade_inercial

    def _calcular_parametro_orbital(self) -> float:
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

# Exemplo de uso
orbita = Orbita(semi_eixo_maior=7000e3, excentricidade=0.001, inclinacao=98.7, raan=120, arg_perigeu=45, anomalia_verdadeira=10)
posicao, velocidade = orbita.calcular_vetor_estado()
print("Posição:", posicao)
print("Velocidade:", velocidade)