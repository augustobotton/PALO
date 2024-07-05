import numpy as np
import matplotlib.pyplot as plt
from src.domain.utilidades_mecanica_orbital.Orbitas.coe_from_sv import coe_from_sv

RA0 = 0.7853868763995575
w0 = 0.5236311392230265
h0 = 55839.27186487775
i0 = 0.48869352872783733,  # Inclinação inicial (rad)
e0 = 0.17136945687836064  # Excentricidade inicial
#TODO arrumar como passar dados iniciais

def elementos_orbitais_resposta(y, t, mu):
    n_times = len(t)

    r = np.zeros(n_times)
    v = np.zeros(n_times)
    h = np.zeros(n_times)
    e = np.zeros(n_times)
    RA = np.zeros(n_times)
    i = np.zeros(n_times)
    w = np.zeros(n_times)
    TA = np.zeros(n_times)

    for j in range(n_times):
        R = y[j, :3]
        V = y[j, 3:]
        r[j] = np.linalg.norm(R)
        v[j] = np.linalg.norm(V)

        _, e[j], i[j], RA[j], w[j], TA[j], h[j] = coe_from_sv(R, V, mu)

    plt.figure(1)

    plt.subplot(2, 1, 1)
    plt.plot(t / 3600, np.rad2deg(RA-RA0))
    plt.title('Variation of Right Ascension')
    plt.xlabel('hours')
    plt.ylabel('${\\it\\Delta\\Omega}$ (deg)')
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.tight_layout()

    plt.subplot(2, 1, 2)
    plt.plot(t / 3600, np.rad2deg(w - w0))
    plt.title('Variation of Argument of Perigee')
    plt.xlabel('hours')
    plt.ylabel('${\\it\\Delta\\omega}$ (deg)')
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.autoscale(enable=True, axis='both', tight=True)
    plt.tight_layout()

    plt.figure(2)

    plt.subplot(3, 1, 1)
    plt.plot(t / 3600, h - h0)
    plt.title('Variation of Angular Momentum')
    plt.xlabel('hours')
    plt.ylabel('${\\it\\Delta h}$ (km$^2$/s)')
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.autoscale(enable=True, axis='both', tight=True)
    plt.tight_layout()

    plt.subplot(3, 1, 2)
    plt.plot(t / 3600, e - e0)
    plt.title('Variation of Eccentricity')
    plt.xlabel('hours')
    plt.ylabel('${\\it\\Delta e}$')
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.autoscale(enable=True, axis='both', tight=True)
    plt.tight_layout()

    plt.subplot(3, 1, 3)
    plt.plot(t / 3600, np.rad2deg(i - i0))
    plt.title('Variation of Inclination')
    plt.xlabel('hours')
    plt.ylabel('${\\it\\Delta i}$ (deg)')
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.autoscale(enable=True, axis='both', tight=True)
    plt.tight_layout()

    plt.show()


