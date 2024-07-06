import pickle

import numpy as np
from scipy.integrate import solve_ivp
from tqdm import tqdm

from src.domain.modelos.foguete.ModeloDinamica import dinamica_foguete
from src.domain.modelos.manobras.parametros_manobra_adquire_gso import ParametrosManobraAdquireOrbitaDeTransferencia
from src.domain.modelos.planeta.Planeta import ModeloPlaneta as Planeta


class Simulacao:
    """
    Classe para simulação de lançamento de foguete com destino a uma órbita específica.

    Atributos:
    - planeta: Objeto Planeta que representa o planeta de lançamento.
    - base_de_lancamento: Objeto que representa a base de lançamento.
    - foguete: Objeto que representa o foguete.
    - condicoes_iniciais: Lista de condições iniciais para a simulação.
    """

    def __init__(self, planeta: Planeta, base_de_lancamento, foguete, condicoes_iniciais):
        """
        Inicializa a simulação com os parâmetros fornecidos.

        :param planeta: Objeto Planeta que representa o planeta de lançamento.
        :param base_de_lancamento: Objeto que representa a base de lançamento.
        :param foguete: Objeto que representa o foguete.
        :param condicoes_iniciais: Lista de condições iniciais para a simulação.
        """
        self.distancia_radial_inicial = foguete.modelo_propulsivo.distancia_radial_inicial
        self.tempo_simulacao = condicoes_iniciais[0]
        self.velocidade_inicial = condicoes_iniciais[1]
        self.angulo_elevacao_inicial = condicoes_iniciais[2]
        self.orbita_alvo = condicoes_iniciais[3]
        self.planeta = planeta
        self.base_de_lancamento = base_de_lancamento
        self.foguete = foguete

        self.inclinacao_orbita = self.orbita_alvo.inclinacao
        self.altitude_alvo = self.orbita_alvo.semi_eixo_maior

        self._configurar_parametros_base()
        self._calcular_condicoes_azimute()

    def _configurar_parametros_base(self):
        """Configura os parâmetros iniciais da base de lançamento."""
        self.altitude_inicial = self.base_de_lancamento.altitude_base_de_lancamento
        self.latitude_inicial = self.base_de_lancamento.latitude_base_de_lancamento
        self.longitude_inicial = self.base_de_lancamento.longitude_base_de_lancamento
        self.comprimento_trilho = self.base_de_lancamento.comprimento_do_trilho

    def _calcular_condicoes_azimute(self):
        """Calcula as condições de azimute iniciais e finais para a simulação."""
        y_t = np.cos(self.inclinacao_orbita) / np.cos(self.latitude_inicial)
        if abs(y_t) > 1:
            print('Não é possível atingir a inclinação a partir da latitude inicial. Calculando a menor possível')
            y_t = np.sign(y_t)

        self.azimute_final = np.arcsin(y_t)
        self.azimute_inicial = self._estimar_azimute_inicial()

        print(f'Condição final de azimute de velocidade inercial (grau): {np.rad2deg(self.azimute_final)}')
        print(f'Condição inicial de azimute de velocidade relativa (grau): {np.rad2deg(self.azimute_inicial)}')

    def _estimar_azimute_inicial(self):
        """
        Estima o azimute inicial baseado nas condições de transferência para a órbita alvo.

        :return: Azimute inicial estimado (rad).
        """
        apogeu_transferencia = self.planeta.raio_equatorial + 250e3
        semi_eixo_maior_transferencia = (self.altitude_alvo + apogeu_transferencia) / 2
        velocidade_transferencia = np.sqrt(
            self.planeta.mut * (2 / apogeu_transferencia - 1 / semi_eixo_maior_transferencia))
        return np.arctan(np.tan(self.azimute_final) - (
                apogeu_transferencia * self.planeta.velocidade_inercial_de_rotacao * np.cos(self.latitude_inicial)) / (
                                 velocidade_transferencia * np.cos(self.azimute_final)))

    def progress_callback(self, t, y, barra):
        """
        Callback para atualização da barra de progresso durante a simulação.

        :param t: Tempo atual da simulação.
        :param y: Estado atual do sistema.
        :param barra: Objeto da barra de progresso.
        :return: Tempo atualizado.
        """
        barra.update(t - barra.n)
        return t

    def simular(self):
        """
        Executa a simulação do lançamento do foguete.

        :return: Vetores de tempo (t) e estados (y) da simulação.
        """
        parametros_apogeu = ParametrosManobraAdquireOrbitaDeTransferencia()

        condicoes_iniciais = [self.velocidade_inicial, self.azimute_inicial, self.angulo_elevacao_inicial,
                              self.distancia_radial_inicial, self.latitude_inicial, self.longitude_inicial]

        opcoes_integracao = {'rtol': 1e-8, 'atol': 1e-10, 'max_step': 1}
        print(f'Simulando por {self.tempo_simulacao} segundos')

        with tqdm(total=self.tempo_simulacao) as bar:
            def fun(t, y):
                bar.update(t - bar.n)
                return dinamica_foguete(t, y, self.base_de_lancamento, self.planeta, self.foguete, parametros_apogeu,
                                        self.orbita_alvo)

            resposta_simulacao = solve_ivp(
                fun,
                (0, self.tempo_simulacao),
                y0=condicoes_iniciais,
                method='RK45',
                t_eval=None,
                **opcoes_integracao
            )

        # Salvar a resposta em um arquivo
        with open('../construtorderesultados/resposta_simulacao.pkl', 'wb') as f:
            pickle.dump(resposta_simulacao, f)

        return resposta_simulacao.t, resposta_simulacao.y
