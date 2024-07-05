import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks

from src.domain.modelos.orbitas.utilidades.vetor_de_estados_para_coe import vetor_de_estados_para_coe


def plota_variacao_elementos_orbitais(y, t, mu, orbita):
    """
    Calcula e plota a variação dos elementos orbitais ao longo do tempo.

    Args:
        y (np.ndarray): Array contendo as posições e velocidades.
        t (np.ndarray): Array de tempos correspondentes às posições e velocidades.
        mu (float): Parâmetro gravitacional padrão do corpo central.
        orbita (Orbita): Instância da classe Orbita contendo os elementos orbitais iniciais.

    Plota:
        - Variação da ascensão reta do nó ascendente.
        - Variação do argumento do periastro.
        - Variação do momento angular.
        - Variação da excentricidade.
        - Variação da inclinação.
    """
    n_times = len(t)

    ra0 = orbita.raan
    w0 = orbita.arg_periastro
    h0 = orbita.quantidade_momento_angular
    i0 = orbita.inclinacao
    e0 = orbita.excentricidade

    r = np.zeros(n_times)
    v = np.zeros(n_times)
    h = np.zeros(n_times)
    e = np.zeros(n_times)
    ra = np.zeros(n_times)
    i = np.zeros(n_times)
    w = np.zeros(n_times)
    ta = np.zeros(n_times)

    for j in range(n_times):
        R = y[j, :3]
        V = y[j, 3:]
        r[j] = np.linalg.norm(R)
        v[j] = np.linalg.norm(V)

        _, e[j], i[j], ra[j], w[j], ta[j], h[j] = vetor_de_estados_para_coe(R, V, mu)

    # Plotando as variações
    plt.figure(1)
    plt.subplot(2, 1, 1)
    plt.plot(t / 3600, np.rad2deg(ra - ra0))
    plt.title('Variação da Ascensão Reta do Nó Ascendente')
    plt.xlabel('horas')
    plt.ylabel('${\\it\\Delta\\Omega}$ (graus)')
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.tight_layout()

    plt.subplot(2, 1, 2)
    plt.plot(t / 3600, np.rad2deg(w - w0))
    plt.title('Variação do Argumento do Periastro')
    plt.xlabel('horas')
    plt.ylabel('${\\it\\Delta\\omega}$ (graus)')
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.tight_layout()

    plt.figure(2)
    plt.subplot(3, 1, 1)
    plt.plot(t / 3600, h - h0)
    plt.title('Variação do Momento Angular')
    plt.xlabel('horas')
    plt.ylabel('${\\it\\Delta h}$ (km$^2$/s)')
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.tight_layout()

    plt.subplot(3, 1, 2)
    plt.plot(t / 3600, e - e0)
    plt.title('Variação da Excentricidade')
    plt.xlabel('horas')
    plt.ylabel('${\\it\\Delta e}$')
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.tight_layout()

    plt.subplot(3, 1, 3)
    plt.plot(t / 3600, np.rad2deg(i - i0))
    plt.title('Variação da Inclinação')
    plt.xlabel('horas')
    plt.ylabel('${\\it\\Delta i}$ (graus)')
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.tight_layout()

    plt.show()


def set_axes_equal(ax):
    """
    Define a mesma escala para os eixos x, y e z.
    """
    limits = np.array([ax.get_xlim3d(), ax.get_ylim3d(), ax.get_zlim3d()])
    centers = np.mean(limits, axis=1)
    max_range = np.max(np.abs(limits[:, 1] - limits[:, 0])) / 2.0

    ax.set_xlim3d([centers[0] - max_range, centers[0] + max_range])
    ax.set_ylim3d([centers[1] - max_range, centers[1] + max_range])
    ax.set_zlim3d([centers[2] - max_range, centers[2] + max_range])


def plot_orbit(y, R, r0):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Criação da esfera representando o planeta
    u, v = np.mgrid[0:2 * np.pi:100j, 0:np.pi:50j]
    x = R * np.cos(u) * np.sin(v)
    y_sphere = R * np.sin(u) * np.sin(v)
    z = R * np.cos(v)
    ax.plot_surface(x, y_sphere, z, color='grey', alpha=0.5)  # Desenha o planeta

    # Plota a trajetória da órbita
    ax.plot(y[:, 0], y[:, 1], y[:, 2], 'k')
    # Linha radial do ponto inicial
    ax.plot([0, r0[0]], [0, r0[1]], [0, r0[2]], 'r')
    # Marca o ponto inicial
    ax.text(y[0, 0], y[0, 1], y[0, 2], 'o', color='red')
    # Marca o ponto final
    ax.text(y[-1, 0], y[-1, 1], y[-1, 2], 'f', color='blue')

    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')
    ax.set_zlabel('Z (km)')
    ax.set_aspect('auto')

    # Ajustar os eixos para terem o mesmo tamanho
    set_axes_equal(ax)

    plt.show()

def constroe_resultados(RE, r0, t, v0, y):
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
    plt.tight_layout()
    plt.show()
    # Impressão de resultados
    print("\nÓrbita Terrestre")
    print("A posição inicial é", r0)
    print("Magnitude =", np.linalg.norm(r0), "km")
    print("A velocidade inicial é", v0)
    print("Magnitude =", np.linalg.norm(v0), "km/s")
    print("Tempo inicial = {:.2f} h. Tempo final = {:.2f} h.".format(t[0] / 3600, t[-1] / 3600))
    print("A altitude mínima é {:.2f} km no tempo = {:.2f} h.".format(rmin - RE, t[imin] / 3600))
    print("A velocidade nesse ponto é {:.2f} km/s.".format(v_at_rmin))
    print("A altitude máxima é {:.2f} km no tempo = {:.2f} h.".format(rmax - RE, t[imax] / 3600))
    print("A velocidade nesse ponto é {:.2f} km/s.".format(v_at_rmax))

