import numpy as np
from scipy.integrate import solve_ivp

from src.domain.modelos.orbitas.propagacao.numerica.dinamicas.dinamica_perturbada_J2 import dinamica_perturbada_J2
from src.domain.modelos.orbitas.utilidades.calculos_orbitais import calcular_periodo_orbital
from src.domain.modelos.orbitas.utilidades.rv_from_r0v0 import rv_from_r0v0


def propagacao_encke(t0, tf, orbita):
    """
    Realiza a propagação da órbita utilizando o método de Encke com a perturbação J2.

    :param t0: Tempo inicial da integração (em segundos).
    :param tf: Tempo final da integração (em segundos).
    :param orbita: Objeto Orbita contendo os parâmetros orbitais iniciais.
    :return: Vetor de tempos (t) e matriz de estados (y) ao longo da propagação.
    """

    R0 = orbita.vetor_posicao
    V0 = orbita.vetor_velocidade
    T0 = calcular_periodo_orbital(orbita.semi_eixo_maior, orbita.mu)
    # Passo de tempo para o procedimento de Encke
    del_t = T0 / 800

    # Início da integração de Encke
    t = t0
    tsave = [t0]
    y = [np.concatenate((R0, V0))]
    del_y0 = np.zeros(6)
    t += del_t



    opts = {'rtol': 1e-8, 'atol': 1e-8}

    # Loop de integração
    while t <= tf + del_t / 2:
        solver = solve_ivp(dinamica_perturbada_J2, [t0, t], del_y0, args=(orbita.mu, R0, V0, t0, 6378.1370, 1082.63e-6),
                           **opts, max_step=del_t)
        z = solver.y.T[-1]
        Rosc, Vosc = rv_from_r0v0(R0, V0, t - t0, orbita.mu)
        R0 = Rosc + z[:3]
        V0 = Vosc + z[3:]
        t0 = t
        tsave.append(t)
        y.append(np.concatenate((R0, V0)))
        t += del_t
        del_y0 = np.zeros(6)
    y = np.array(y)

    return tsave, y
