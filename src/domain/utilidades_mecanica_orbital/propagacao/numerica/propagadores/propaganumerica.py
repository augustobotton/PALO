import numpy as np
from scipy.integrate import solve_ivp

from src.domain.utilidades_mecanica_orbital.Utilidades.plotaOrbita import plot_orbit
from src.domain.utilidades_mecanica_orbital.propagacao.numerica.dinamicas.dinamicas_orbitais import \
    dinamica_lambert


def propagacao_numerica(ti, tf, r0, v0):
    hours = 3600
    R = 6378  # Planet radius (km)
    y0 = np.array(r0 + v0)
    opts = {'rtol': 1e-8, 'atol': 1e-8}

    args = [398600.4418]
    sol = solve_ivp(dinamica_lambert, [ti, tf], y0, args=args, **opts)

    t = sol.t
    y = sol.y.T

    r = np.linalg.norm(y[:, :3], axis=1)
    rmax = np.max(r)
    rmin = np.min(r)
    imax = np.argmax(r)
    imin = np.argmin(r)

    v_at_rmax = np.linalg.norm(y[imax, 3:6])
    v_at_rmin = np.linalg.norm(y[imin, 3:6])

    print("\nEarth Orbit")
    print("The initial position is", r0)
    print("Magnitude =", np.linalg.norm(r0), "km")
    print("The initial velocity is", v0)
    print("Magnitude =", np.linalg.norm(v0), "km/s")
    print("Initial time = {:.2f} h. Final time = {:.2f} h.".format(ti / hours, tf / hours))
    print("The minimum altitude is {:.2f} km at time = {:.2f} h.".format(rmin - R, t[imin] / hours))
    print("The speed at that point is {:.2f} km/s.".format(v_at_rmin))
    print("The maximum altitude is {:.2f} km at time = {:.2f} h.".format(rmax - R, t[imax] / hours))
    print("The speed at that point is {:.2f} km/s.".format(v_at_rmax))
    plot_orbit(y, R, r0)


