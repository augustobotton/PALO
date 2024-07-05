import numpy as np
from scipy.integrate import solve_ivp

from src.domain.utilidades_mecanica_orbital.Orbitas.rv_from_coe import sv_from_coe
from src.domain.utilidades_mecanica_orbital.Utilidades.calculos_orbitais import calcular_periodo_orbital
from src.domain.utilidades_mecanica_orbital.propagacao.numerica.dinamicas.dinamicas_orbitais import dinamica_perturbada_J2
from src.domain.utilidades_mecanica_orbital.Utilidades.rv_from_r0v0 import rv_from_r0v0




def propagacao_encke(t0, tf, orbita):

    coe = orbita.retorna_parametros()

    R0, V0 = sv_from_coe(coe, orbita.mu)

    T0 = calcular_periodo_orbital(orbita.semi_eixo_maior, orbita.mu)

    # Time step for Encke procedure
    del_t = T0 / 100

    # Begin the Encke integration
    t = t0
    tsave = [t0]
    y = np.hstack((R0, V0))
    y0 = np.zeros((6))
    t += del_t

    args = [orbita.mu, R0, V0, t0, 6378.1370, 0.001082630]

    opts = {'rtol': 1e-8, 'atol': 1e-8, 'max_step': del_t}
    print('ROVO antes de integrar',R0,V0)
    # Integration loop
    while t <= tf + del_t / 2:

        solver = solve_ivp(dinamica_perturbada_J2, [t0, t], y0, args=args, **opts)

        # Compute the osculating state vector at time t
        Rosc, Vosc = rv_from_r0v0(R0, V0, t - t0, orbita.mu)

        # Rectify
        R0 = Rosc + solver.y[:3, -1]
        V0 = Vosc + solver.y[3:, -1]
        t0 = t
        # Prepare for next time step
        tsave.append(t)
        t += del_t
        y = np.vstack((y, np.hstack((R0, V0))))
        y0 = np.zeros((6))

    t = np.array(tsave)
    return t, y
