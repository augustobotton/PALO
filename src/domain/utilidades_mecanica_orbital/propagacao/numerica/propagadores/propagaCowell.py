import numpy as np
from scipy.integrate import solve_ivp
from scipy.signal import find_peaks
from matplotlib import pyplot as plt
from src.domain.utilidades_mecanica_orbital import Utilidades
from src.domain.utilidades_mecanica_orbital.Orbitas.ModeloOrbita import Orbita
from src.domain.utilidades_mecanica_orbital.Orbitas.rv_from_coe import sv_from_coe
from src.domain.utilidades_mecanica_orbital.Utilidades.plotaOrbita import plot_orbit
from src.domain.utilidades_mecanica_orbital.propagacao.numerica.dinamicas.dinamica_p_arrasto import \
    dinamica_perturbada_arrasto

def terminate(t, y, *args):
    r = np.linalg.norm(y[:3])
    altitude = r - 6378
    return altitude - 100

terminate.terminal = True
terminate.direction = -1  # O evento é acionado quando a função diminui e passa por zero

def propaga_cowell(ti, tf, orbita: Orbita, atributos_dinamica):
    hours = 3600
    coe = orbita.retorna_parametros()
    r0, v0 = sv_from_coe(coe, orbita.mu)
    y0 = np.concatenate((np.array(r0), np.array(v0)))
    T = Utilidades.calculos_orbitais.calcular_periodo_orbital(coe[0], orbita.mu)
    opts = {'rtol': 1e-8, 'atol': 1e-8, 'first_step': T/10000}
    nout = 40000  # Número de pontos de saída
    t_eval = np.linspace(ti, tf, nout)  # Intervalo de tempo de integração
    sol = solve_ivp(dinamica_perturbada_arrasto, [ti, tf], y0, events=terminate, args=atributos_dinamica, **opts)
    terra = atributos_dinamica[0]
    RE = terra.raio_equatorial
    t = sol.t
    y = sol.y.T

    r = np.linalg.norm(y[:, :3], axis=1)
    v = np.linalg.norm(y[:, 3:6], axis=1)
    # Máximo e mínimo globais da altitude
    rmax = np.max(r)
    rmin = np.min(r)
    imax = np.argmax(r)
    imin = np.argmin(r)

    # Velocidade nos pontos de máxima e mínima altitude
    v_at_rmax = v[imax]
    v_at_rmin = v[imin]

    # Altitude a cada instante
    altitude = r - RE

    # Encontrar os máximos e mínimos locais de altitude
    imax_local, _ = find_peaks(altitude)
    imin_local, _ = find_peaks(-altitude)
    maxima = np.column_stack((t[imax_local], altitude[imax_local]))
    minima = np.column_stack((t[imin_local], -altitude[imin_local]))

    # Ordenar por tempo
    apogee = maxima[maxima[:, 0].argsort()]
    perigee = minima[minima[:, 0].argsort()]

    # Plotagem dos máximos e mínimos locais
    plt.figure(1)
    plt.plot(apogee[:, 0], apogee[:, 1], 'b', linewidth=2, label='Apogee')
    plt.plot(perigee[:, 0], perigee[:, 1], 'r', linewidth=2, label='Perigee')
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.xlabel('Time (seconds)')
    plt.ylabel('Altitude (km)')
    plt.ylim([0, 1000])
    plt.legend()

    # Impressão de resultados
    print("\nEarth Orbit")
    print("The initial position is", r0)
    print("Magnitude =", np.linalg.norm(r0), "km")
    print("The initial velocity is", v0)
    print("Magnitude =", np.linalg.norm(v0), "km/s")
    print("Initial time = {:.2f} h. Final time = {:.2f} h.".format(t[0] / 3600, t[-1] / 3600))
    print("The minimum altitude is {:.2f} km at time = {:.2f} h.".format(rmin - RE, t[imin] / 3600))
    print("The speed at that point is {:.2f} km/s.".format(v_at_rmin))
    print("The maximum altitude is {:.2f} km at time = {:.2f} h.".format(rmax - RE, t[imax] / 3600))
    print("The speed at that point is {:.2f} km/s.".format(v_at_rmax))

    return t, y

