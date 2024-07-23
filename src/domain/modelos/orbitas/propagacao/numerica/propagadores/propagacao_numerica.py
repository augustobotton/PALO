import numpy as np
from scipy.integrate import solve_ivp

from src.domain.modelos.orbitas.propagacao.numerica.dinamicas.dinamica_lambert import \
    dinamica_lambert


def propagacao_numerica(ti, tf, orbita):

    r0, v0 = orbita.vetor_posicao, orbita.vetor_velocidade
    y0 = np.concatenate((r0, v0))
    opts = {'rtol': 1e-8, 'atol': 1e-8}

    args = [398600.4418]
    sol = solve_ivp(dinamica_lambert, [ti, tf], y0, args=args, **opts)

    t = sol.t
    y = sol.y.T

    return t, y


