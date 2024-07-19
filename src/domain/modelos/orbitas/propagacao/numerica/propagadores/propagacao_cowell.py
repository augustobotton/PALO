import numpy as np
from scipy.integrate import solve_ivp

from src.domain.modelos.orbitas.Orbita import Orbita
from src.domain.modelos.orbitas.propagacao.numerica.dinamicas.dinamica_perturbada_arrasto import \
    dinamica_perturbada_arrasto
from src.domain.modelos.orbitas.utilidades.calculos_orbitais import calcular_periodo_orbital
from src.domain.modelos.orbitas.utilidades.coe_para_vetor_de_estados import coe_para_vetor_de_estados


def terminate(t, y, *args):
    """
    Função de término para a integração numérica, baseada na altitude da órbita.

    Parameters:
    t : float
        Tempo atual da integração.
    y : array_like
        Vetor de estado atual.

    Returns:
    float
        A diferença entre a altitude atual e 100 km.
    """
    altitude_limite = args[4]
    r = np.linalg.norm(y[:3])
    altitude = r - 6378
    return altitude - altitude_limite


terminate.terminal = True
terminate.direction = -1  # O evento é acionado quando a função diminui e passa por zero


def propaga_cowell(ti, tf, orbita: Orbita, atributos_dinamica):
    """
    Propaga a órbita utilizando o método de Cowell com perturbações de arrasto.

    Parameters:
    ti : float
        Tempo inicial da integração (em segundos).
    tf : float
        Tempo final da integração (em segundos).
    orbita : Orbita
        Objeto Orbita contendo os parâmetros orbitais iniciais.
    atributos_dinamica : tuple
        Atributos necessários para a dinâmica perturbada (ex. modelo da Terra).

    Returns:
    t : ndarray
        Vetor de tempos da solução.
    y : ndarray
        Matriz de estados solucionados em cada instante de tempo.
    """

    coe = orbita.retorna_parametros()
    r0, v0 = coe_para_vetor_de_estados(coe, orbita.mu)
    y0 = np.concatenate((np.array(r0), np.array(v0)))
    T = calcular_periodo_orbital(coe[0], orbita.mu)
    opts = {'rtol': 1e-8, 'atol': 1e-8, 'first_step': T / 10000}
    nout = 40000  # Número de pontos de saída
    t_eval = np.linspace(ti, tf, nout)  # Intervalo de tempo de integração
    sol = solve_ivp(dinamica_perturbada_arrasto, [ti, tf], y0, events=terminate, args=atributos_dinamica, **opts)
    t = sol.t
    y = sol.y.T

    return t, y

