import pickle
from tqdm import tqdm
import numpy as np
from scipy.integrate import solve_ivp

from src.domain.modelos.foguete.ModeloDinamica import dinamica_foguete
from src.domain.modelos.manobras.parametros_manobra_adquire_gso import ParametrosManobraAdquireOrbitaDeTransferencia


class Simulacao:
    def __init__(self, planeta, base_de_lancamento, foguete, condicoes_iniciais):
        self.tempo_simulacao = condicoes_iniciais[0]
        self.velocidade_inicial = condicoes_iniciais[1]
        self.angulo_elevacao_inicial = np.radians(condicoes_iniciais[2])
        self.orbita_alvo = condicoes_iniciais[3]
        self.planeta = planeta
        self.base_de_lancamento = base_de_lancamento
        self.foguete = foguete

        self.inclinacao_orbita = np.radians(self.orbita_alvo.inclinacao)
        self.altitude_alvo = self.orbita_alvo.semi_eixo_maior

        self._configurar_parametros_planeta()
        self._configurar_parametros_base()
        self._configurar_parametros_orbita()
        self._calcular_condicoes_azimute()

    def _configurar_parametros_planeta(self):
        self.raio_equatorial = self.planeta.raio_equatorial
        self.velocidade_rotacao = self.planeta.velocidade_inercial_de_rotacao
        self.gravidade = self.planeta.gravidade
        self.constante_gravitacional = self.planeta.mut
        self.j2 = self.planeta.J2
        self.j3 = self.planeta.J3
        self.j4 = self.planeta.J4

    def _configurar_parametros_base(self):
        self.altitude_inicial = self.base_de_lancamento.altitude_base_de_lancamento
        self.latitude_inicial = self.base_de_lancamento.latitude_base_de_lancamento
        self.longitude_inicial = self.base_de_lancamento.longitude_base_de_lancamento
        self.comprimento_trilho = self.base_de_lancamento.comprimento_do_trilho

    def _configurar_parametros_orbita(self):
        self.velocidade_geo_sincrona = np.sqrt(self.constante_gravitacional / self.altitude_alvo)
        self.distancia_radial_inicial = self.raio_equatorial + self.altitude_inicial

    def _calcular_condicoes_azimute(self):
        y_t = np.cos(self.inclinacao_orbita) / np.cos(self.latitude_inicial)
        if abs(y_t) > 1:
            print('Nao eh possivel atingir a inclinacao a partir da latitude inicial. Calculando a menor possivel')
            y_t = np.sign(y_t)

        self.azimute_final = np.arcsin(y_t)
        self.azimute_inicial = self._estimar_azimute_inicial()

        print(f'Condicao final de azimute de velocidade inercial (grau): {np.degrees(self.azimute_final)}')
        print(f'Condicao inicial de azimute de velocidade relativa (grau): {np.degrees(self.azimute_inicial)}')

    def _estimar_azimute_inicial(self):
        apogeu_transferencia = self.raio_equatorial + 250e3
        semi_eixo_maior_transferencia = (self.altitude_alvo + apogeu_transferencia) / 2
        velocidade_transferencia = np.sqrt(
            self.constante_gravitacional * (2 / apogeu_transferencia - 1 / semi_eixo_maior_transferencia))
        return np.arctan(np.tan(self.azimute_final) - (
                apogeu_transferencia * self.velocidade_rotacao * np.cos(self.latitude_inicial)) / (
                                 velocidade_transferencia * np.cos(self.azimute_final)))

    def progress_callback(self, t, y, bar):
        bar.update(t - bar.n)
        return t

    def simular(self):
        parametros_apogeu = ParametrosManobraAdquireOrbitaDeTransferencia()

        condicoes_iniciais = [self.velocidade_inicial, self.azimute_inicial, self.angulo_elevacao_inicial,
                              self.distancia_radial_inicial, self.latitude_inicial, self.longitude_inicial]

        opcoes_integracao = {'rtol': 1e-8, 'atol': 1e-10, 'max_step': 1}
        print(f'Simulando por {self.tempo_simulacao} segundos')
        # Create a progress bar
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

        # Save the response to a file
        with open('resposta_simulacao.pkl', 'wb') as f:
            pickle.dump(resposta_simulacao, f)

        return resposta_simulacao.t, resposta_simulacao.y


