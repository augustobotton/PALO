import numpy as np
from scipy.interpolate import pchip_interpolate


#TODO: OK!
class ModeloAerodinamico:
    def __init__(self):
        self.altitude = None
        self.numero_de_knudsen = None
        self.numero_de_mach = None
        self.temperatura = None
        self.constante_do_gas_ideal = None
        self.velocidade = None

    def calcula(self) -> float:
        """
        Calcula o coeficiente de arrasto com base nos parâmetros fornecidos.

        :return: Coeficiente de arrasto
        """
        coeficiente_de_arrasto_em_escoamento_livre_molecular = self._computa_coef_de_arrasto_em_escoamento_livre_molecular()
        coeficiente_de_arrasto_interpolado = self._interpola_coeficiente_de_arrasto(self.numero_de_mach)
        if self.altitude < 2000e3:
            if self.numero_de_knudsen < 0.0146:
                coeficiente_de_arrasto = coeficiente_de_arrasto_interpolado
            elif self.numero_de_knudsen < 14.5:
                coeficiente_de_arrasto = coeficiente_de_arrasto_interpolado + (
                        coeficiente_de_arrasto_em_escoamento_livre_molecular - coeficiente_de_arrasto_interpolado) * (
                                                 (1 / 3) * np.log10(
                                             self.numero_de_knudsen / np.sin(30 * np.pi / 180)) + 0.05113)
            else:
                coeficiente_de_arrasto = coeficiente_de_arrasto_em_escoamento_livre_molecular
        else:
            coeficiente_de_arrasto = 0

        return coeficiente_de_arrasto

    def _computa_coef_de_arrasto_em_escoamento_livre_molecular(self) -> float:
        """
        Computa o coeficiente de arrasto em escoamento livre molecular.

        :return: Coeficiente de arrasto em escoamento livre molecular
        """
        return 1.75 + np.sqrt(np.pi) / (
                2 * (self.velocidade / np.sqrt(2 * self.constante_do_gas_ideal * self.temperatura)))

    @staticmethod
    def _interpola_coeficiente_de_arrasto(numero_de_mach: float) -> float:
        """
        Interpola o coeficiente de arrasto com base no número de Mach.

        :param numero_de_mach: Número de Mach
        :return: Coeficiente de arrasto interpolado
        """
        coeficiente_de_arrasto_por_mach = np.array([[0.0000, 0.4736],
                                                    [0.1609, 0.4736],
                                                    [0.3448, 0.4764],
                                                    [0.4828, 0.4875],
                                                    [0.5261, 0.5001],
                                                    [0.6424, 0.5141],
                                                    [0.7091, 0.5389],
                                                    [0.7521, 0.5646],
                                                    [0.8463, 0.5905],
                                                    [0.8736, 0.6153],
                                                    [0.9065, 0.6667],
                                                    [0.9378, 0.6375],
                                                    [0.9492, 0.7029],
                                                    [0.9530, 0.7217],
                                                    [0.9721, 0.7858],
                                                    [0.9759, 0.8046],
                                                    [0.9949, 0.8680],
                                                    [0.9987, 0.8866],
                                                    [1.0049, 1.0313],
                                                    [1.0115, 0.9833],
                                                    [1.0178, 0.9502],
                                                    [1.0247, 1.0133],
                                                    [1.0345, 0.9667],
                                                    [1.1264, 1.0125],
                                                    [1.1494, 0.9958],
                                                    [1.2184, 0.9833],
                                                    [1.2644, 0.9653],
                                                    [1.2656, 0.9569],
                                                    [1.3037, 0.9398],
                                                    [1.3494, 0.9193],
                                                    [1.3952, 0.8998],
                                                    [1.4409, 0.8797],
                                                    [1.4866, 0.8597],
                                                    [1.5324, 0.8415],
                                                    [1.5862, 0.8194],
                                                    [1.6696, 0.7925],
                                                    [1.7382, 0.7706],
                                                    [1.8068, 0.7491],
                                                    [1.8754, 0.7277],
                                                    [1.9770, 0.6944],
                                                    [2.1842, 0.6597],
                                                    [2.2757, 0.6445],
                                                    [2.3557, 0.6315],
                                                    [2.6073, 0.5941],
                                                    [2.7331, 0.5785],
                                                    [2.8703, 0.5596],
                                                    [2.9885, 0.5458],
                                                    [3.1695, 0.5396],
                                                    [3.4306, 0.5272],
                                                    [3.6822, 0.5171],
                                                    [3.9719, 0.5109],
                                                    [4.1853, 0.5060],
                                                    [4.4369, 0.5019],
                                                    [4.6885, 0.4979],
                                                    [4.9545, 0.4920],
                                                    [5.1821, 0.4923],
                                                    [5.4432, 0.4921],
                                                    [5.6948, 0.4921],
                                                    [5.9463, 0.4908],
                                                    [6.1979, 0.4908],
                                                    [6.4495, 0.4901],
                                                    [6.7011, 0.4894],
                                                    [6.9526, 0.4893],
                                                    [7.2042, 0.4881],
                                                    [7.4558, 0.4880],
                                                    [7.7073, 0.4872],
                                                    [7.9589, 0.4866],
                                                    [8.2105, 0.4865],
                                                    [8.4621, 0.4853],
                                                    [8.7136, 0.4852],
                                                    [8.9652, 0.4845],
                                                    [9.2168, 0.4839],
                                                    [9.4683, 0.4837],
                                                    [9.7199, 0.4825],
                                                    [10.2231, 0.4823],
                                                    [10.4746, 0.4825],
                                                    [10.7262, 0.4825],
                                                    [10.9778, 0.4825],
                                                    [11.2294, 0.4825],
                                                    [11.4809, 0.4825],
                                                    [11.7325, 0.4824]])
        return pchip_interpolate(coeficiente_de_arrasto_por_mach[:, 0], coeficiente_de_arrasto_por_mach[:, 1],
                                 numero_de_mach)

    def aerodinamica_multiplos_estagios(self, tempo: float, velocidade: float, area_de_referencia: list,
                                        tempo_limite_separacao: list,
                                        densidade_do_ar: float):
        """
        Calcula as forças aerodinâmicas para foguetes de múltiplos estágios.

        :param tempo: Tempo atual
        :param velocidade: Velocidade do foguete
        :param densidade_do_ar: Densidade do ar
        :return: Tupla contendo o arrasto, a força de sustentação lateral e a força de sustentação
        """

        fator_correcao_arrasto = 1.28

        # Calcular o coeficiente de arrasto
        coeficiente_arrasto = self.calcula()
        coeficiente_arrasto_ajustado = fator_correcao_arrasto * coeficiente_arrasto

        # Determinar a área de referência com base no estágio
        area_do_estagio = self._define_area_referencia(tempo, tempo_limite_separacao, area_de_referencia)

        # Calcular as forças
        arrasto = 0.5 * densidade_do_ar * velocidade ** 2 * area_do_estagio * coeficiente_arrasto_ajustado
        forca_sustentacao_lateral = 0
        forca_sustentacao = 0

        return arrasto, forca_sustentacao_lateral, forca_sustentacao

    def atualizar_parametros(self, **kwargs):
        for key, value in kwargs.items():
            if hasattr(self, key):
                # Tente converter para float se for um dos parâmetros numéricos
                if key in ['velocidade', 'constante_do_gas_ideal', 'temperatura']:
                    try:
                        value = float(value)
                    except ValueError:
                        raise ValueError(f"Valor inválido para {key}: {value} não pode ser convertido para float.")
                setattr(self, key, value)
            else:
                raise ValueError(f"Parâmetro '{key}' não é um atributo válido do ModeloAerodinamico")

    def _define_area_referencia(self, tempo: float, tempo_limite_separacao: list, area_de_referencia: list) -> float:
        """
        Calcula a área de referência com base no estágio do foguete.

        :param tempo: Tempo atual
        :param tempo_limite_separacao: Lista de tempos limite de separação
        :param area_de_referencia: Lista de áreas de referência
        :return: Área de referência para o estágio atual
        """
        n_estagios = len(tempo_limite_separacao)
        for i in range(n_estagios):
            if tempo <= tempo_limite_separacao[i]:
                return area_de_referencia[i]
        return area_de_referencia[-1]  # Padrão para a área do último estágio
