import numpy as np

from src.domain.modelos.foguete.propulsao.ModeloPropulsivo import ModeloPropulsivo
from src.domain.modelos.orbitas.Orbita import Orbita
from src.domain.modelos.orbitas.utilidades import calculos_orbitais
from src.domain.modelos.orbitas.utilidades.funcoes_conversao import componentes_vel_relativa_para_inercial


class ParametrosManobraAdquireOrbitaDeTransferencia():

    def __init__(self):
        self.achou_apogeu = 0
        self.sinal_phi_inercial = 0

    def parametros_manobra_adquire_gso(self, tempo, m, X, orbita_transferencia: Orbita,
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

        # Desmembrar o vetor de estado
        V = X[0]
        A = X[1]
        phi = X[2]
        r = X[3]
        delta = X[4]

        # Vetor velocidade inercial
        Vi, phii, _ = componentes_vel_relativa_para_inercial(V, phi, A, we, r, delta)

        # Realização de uma sequência de testes para verificar a ocorrência do apogeu da órbita GTO.
        # Quando ele ocorre, determina os parâmetros da manobra.
        if r > 0.9 * orbita_transferencia.semi_eixo_maior:
            if not self.achou_apogeu:
                if np.sign(phii) != self.sinal_phi_inercial:
                    # Se o sinal for diferente, phii passou por zero, o foguete chegou no apogeu
                    self.achou_apogeu = True  # Encontrou o apogeu, ao setar essa variável, só vai entrar aqui uma vez
                    print('manobra')
                    modelo_propulsivo.tempos_de_ignicao[3] = tempo + 1  # Guarda o tempo de ignição
                    dv_transferencia = velocidade_orbital_desejada - Vi
                    mp32 = (m * np.exp(dv_transferencia/ (modelo_propulsivo.impulso_especifico[2] * g)) - m) / np.exp(
                        dv_transferencia / (modelo_propulsivo.impulso_especifico[2] * g))

                    # Duração da queima necessária
                    Tq32 = modelo_propulsivo.duracao_total_de_queima_do_terceiro_estagio * mp32 / \
                           modelo_propulsivo.massa_propelente_estagios[3]

                    if modelo_propulsivo.duracao_da_primeira_queima_3_estagio + Tq32 > modelo_propulsivo.duracao_total_de_queima_do_terceiro_estagio:
                        Tq32 = modelo_propulsivo.duracao_total_de_queima_do_terceiro_estagio - \
                               modelo_propulsivo.duracao_da_primeira_queima_3_estagio

                    modelo_propulsivo.tempos_de_fim_de_queima[3] = modelo_propulsivo.tempos_de_ignicao[3] + Tq32
                    modelo_propulsivo.tempos_de_separacao[2] = modelo_propulsivo.tempos_de_fim_de_queima[3] + \
                                                               modelo_propulsivo.tempos_de_separacao[2]

                    print('tempos de ignicao', modelo_propulsivo.tempos_de_ignicao)
                    print('tempos de fim de queima', modelo_propulsivo.tempos_de_fim_de_queima)
                    print('tempos de separacao', modelo_propulsivo.tempos_de_separacao)


        self.sinal_phi_inercial = np.sign(
            phii)  # Guarda o sinal de phi inercial para verificar mudança na próxima iteração


"ti  = lista tempo de ignicao"
"tq = lista tempo de fim de queima"
"ts = lista tempo de separação"
