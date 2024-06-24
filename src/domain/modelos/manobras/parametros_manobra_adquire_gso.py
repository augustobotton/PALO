import numpy as np

from src.domain.utilidades_mecanica_orbital.orbitalUtils import calculos_orbitais
from src.domain.utilidades_mecanica_orbital.Orbitas import ModeloOrbita
from src.domain.utilidades_mecanica_orbital.orbitalUtils.Converte import Vrel2Vine
from src.domain.modelos.foguete.ModeloPropulsivo import ModeloPropulsivo


class ParametrosManobraAdquireOrbitaDeTransferencia():

    def __init__(self):
        self.achou_apogeu = 0
        self.sinal_phi_inercial = np.sign

    def parametros_manobra_adquire_gso(self, t, m, X, orbita_transferencia: ModeloOrbita,
                                       modelo_propulsivo: ModeloPropulsivo, planeta):
        """
        Função para calcular os parâmetros da manobra mono impulsiva de aquisição de órbita GSO.
        Deve ser chamada no final da função de dinâmica, pois precisa rodar ao longo do tempo da simulação,
        verificando quando ocorre o apogeu da órbita GTO e determinando os parâmetros para imprimir o impulso
        de velocidade de circularização. Esta função não retorna saídas, mas atualiza parâmetros por variáveis globais,
        que serão usados na função de propulsão.

        Parâmetros:
        - t (s): tempo
        - m (kg): massa do foguete
        - X: vetor de estado [V, A, phi, r, delta]

        """

        we = planeta.velocidade_inercial_de_rotacao
        g = planeta.gravidade

        velocidade_orbital_desejada = calculos_orbitais.calcula_velocidade_orbital(planeta.mut,
                                                                                   orbita_transferencia.semi_eixo_maior)

        # Desmenbra o vetor de estado
        V = X[0]
        A = X[1]
        phi = X[2]
        r = X[3]
        delta = X[4]

        # Vetor velocidade inercial
        Vi, phii, _ = Vrel2Vine(V, phi, A, we, r, delta)

        # Realização de uma sequência de testes para verificar a ocorrência do apogeu da órbita GTO.
        # Quando ele ocorre, determina os parâmetros da manobra.
        if r > 0.9 * orbita_transferencia.semi_eixo_maior:
            if not self.achou_apogeu:
                if np.sign(phii) != self.sinal_phi_inercial:
                    # Se o sinal for diferente, phii passou por zero, o foguete chegou no apogeu
                    self.achou_apogeu = True  # Encontrou o apogeu, ao setar essa variável, só vai entrar aqui uma vez
                    modelo_propulsivo.tempos_de_ignicao[-1] = t  # Guarda o tempo de ignição
                    dv_transferencia = velocidade_orbital_desejada - Vi
                    mp32 = (m * np.exp(dv_transferencia / (modelo_propulsivo.impulso_especifico[2] * g)) - m) / np.exp(
                        dv_transferencia / (modelo_propulsivo.impulso_especifico[2] * g))
                    modelo_propulsivo.tempos_de_fim_de_queima[-1] = modelo_propulsivo.tempos_de_ignicao[
                                                                        -1] + mp32 / modelo_propulsivo.massa_total_propelente_terceiro_estagio
                    if modelo_propulsivo.tempos_de_fim_de_queima[2] + modelo_propulsivo.tempos_de_fim_de_queima[
                        3] > modelo_propulsivo.duracao_total_de_queima_do_terceiro_estagio:
                        modelo_propulsivo.tempos_de_fim_de_queima[
                            3] = modelo_propulsivo.duracao_total_de_queima_do_terceiro_estagio - \
                                 modelo_propulsivo.tempos_de_fim_de_queima[2]

                    modelo_propulsivo.tempos_de_fim_de_queima[-1] = modelo_propulsivo.tempos_de_fim_de_queima[-1] + \
                                                                    modelo_propulsivo.tempos_de_fim_de_queima[
                                                                        3]
                    modelo_propulsivo.tempos_de_separacao[-1] = modelo_propulsivo.tempos_de_fim_de_queima[-1] + \
                                                                modelo_propulsivo.tempos_de_separacao[2]
                    modelo_propulsivo.calcular_tempos()
        self.sinal_phi_inercial = np.sign(
            phii)  # Guarda o sinal de phi inercial para verificar mudança na próxima iteração
        return 0


"ti  = lista tempo de ignicao"
"tq = lista tempo de fim de queima"
"ts = lista tempo de separação"
