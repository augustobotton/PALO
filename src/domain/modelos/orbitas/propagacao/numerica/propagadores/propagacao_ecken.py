import numpy as np
from scipy.integrate import solve_ivp

from src.domain.modelos.orbitas.propagacao.numerica.dinamicas.dinamica_perturbada_J2 import dinamica_perturbada_J2
from src.domain.modelos.orbitas.utilidades.calculos_orbitais import calcular_periodo_orbital
from src.domain.modelos.orbitas.utilidades.coe_para_vetor_de_estados import coe_para_vetor_de_estados
from src.domain.modelos.orbitas.utilidades.rv_from_r0v0 import rv_from_r0v0


def propagacao_encke(t0, tf, orbita):
    """
    Realiza a propagação da órbita utilizando o método de Encke com a perturbação J2.

    :param t0: Tempo inicial da integração (em segundos).
    :param tf: Tempo final da integração (em segundos).
    :param orbita: Objeto Orbita contendo os parâmetros orbitais iniciais.
    :return: Vetor de tempos (t) e matriz de estados (y) ao longo da propagação.
    """
    coe = orbita.retorna_parametros()

    R0, V0 = coe_para_vetor_de_estados(coe, orbita.mu)
    T0 = calcular_periodo_orbital(orbita.semi_eixo_maior, orbita.mu)

    # Passo de tempo para o procedimento de Encke
    del_t = T0 / 100

    # Início da integração de Encke
    t = t0
    tsave = [t0]
    y = np.hstack((R0, V0))
    y0 = np.zeros((6))
    t += del_t

    args = [orbita.mu, R0, V0, t0, 6378.1370, 0.001082630]

    opts = {'rtol': 1e-8, 'atol': 1e-8, 'max_step': del_t}


    # Loop de integração
    while t <= tf + del_t / 2:
        solver = solve_ivp(dinamica_perturbada_J2, [t0, t], y0, args=args, **opts)

        # Computar o vetor de estado osculador no tempo t
        Rosc, Vosc = rv_from_r0v0(R0, V0, t - t0, orbita.mu)

        # Retificar
        R0 = Rosc + solver.y[:3, -1]
        V0 = Vosc + solver.y[3:, -1]
        t0 = t

        # Preparar para o próximo passo de tempo
        tsave.append(t)
        t += del_t
        y = np.vstack((y, np.hstack((R0, V0))))
        y0 = np.zeros(6)

    t = np.array(tsave)
    return t, y
