import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks

from src.domain.modelos.orbitas.utilidades.vetor_de_estados_para_coe import vetor_de_estados_para_coe


def plota_variacao_elementos_orbitais(tsave, y, orbita):
    """
    Calcula e plota a variação dos elementos orbitais ao longo do tempo.

    Args:
        y (np.ndarray): Array contendo as posições e velocidades.
        t (np.ndarray): Array de tempos correspondentes às posições e velocidades.
        orbita (Orbita): Instância da classe Orbita contendo os elementos orbitais iniciais.

    Plota:
        - Variação da ascensão reta do nó ascendente.
        - Variação do argumento do periastro.
        - Variação do momento angular.
        - Variação da excentricidade.
        - Variação da inclinação.
    """

    RA0 = orbita.raan
    w0 = orbita.arg_periastro
    h0 = orbita.quantidade_momento_angular
    e0 = orbita.excentricidade
    i0 = orbita.inclinacao

    # Extração dos elementos orbitais
    n_times = len(tsave)
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
        coe = vetor_de_estados_para_coe(R, V, 398600)

        e[j] = coe[1]
        i[j] = coe[2]
        RA[j] = coe[3]
        w[j] = coe[4]
        TA[j] = coe[5]
        h[j] = coe[6]

    plt.figure(figsize=(10, 8))

    # Variação da Ascensão Reta
    plt.subplot(2, 1, 1)
    plt.plot(np.array(tsave) / 3600, np.rad2deg(RA - RA0))
    plt.title('Variação da Ascensão Reta')
    plt.xlabel('Horas')
    plt.ylabel(r'$\Delta\Omega$ (graus)')
    plt.grid(True)

    # Variação do Argumento do Perigeu
    plt.subplot(2, 1, 2)
    plt.plot(np.array(tsave) / 3600, np.rad2deg(w - w0))
    plt.title('Variação do Argumento do Perigeu')
    plt.xlabel('Horas')
    plt.ylabel(r'$\Delta\omega$ (graus)')
    plt.grid(True)

    plt.tight_layout()
    plt.show()

    plt.figure(figsize=(10, 12))

    # Variação do Momento Angular
    plt.subplot(3, 1, 1)
    plt.plot(np.array(tsave) / 3600, h - h0)
    plt.title('Variação do Momento Angular')
    plt.xlabel('Horas')
    plt.ylabel(r'$\Delta h$ (km²/s)')
    plt.grid(True)

    # Variação da Excentricidade
    plt.subplot(3, 1, 2)
    plt.plot(np.array(tsave) / 3600, e - e0)
    plt.title('Variação da Excentricidade')
    plt.xlabel('Horas')
    plt.ylabel(r'$\Delta e$')
    plt.grid(True)

    # Variação da Inclinação
    plt.subplot(3, 1, 3)
    plt.plot(np.array(tsave) / 3600, np.rad2deg(i - i0))
    plt.title('Variação da Inclinação')
    plt.xlabel('Horas')
    plt.ylabel(r'$\Delta i$ (graus)')
    plt.grid(True)

    plt.tight_layout()
    plt.show()


def set_axes_equal(ax):

    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)

    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.5 * max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])

def plota_orbita(y, raio_equatorial, r0):

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    # Criação da esfera representando o planeta
    u, v = np.mgrid[0:2 * np.pi:100j, 0:np.pi:50j]
    x = raio_equatorial * np.cos(u) * np.sin(v)
    y_sphere = raio_equatorial * np.sin(u) * np.sin(v)
    z = raio_equatorial * np.cos(v)
    ax.plot_surface(x, y_sphere, z, color='lightgrey', alpha=0.6, shade=True)  # Desenha o planeta

    # Plota a trajetória da órbita
    ax.plot(y[:, 0], y[:, 1], y[:, 2], color='#1f77b4', linewidth=0.5, label='Órbita')  # Cor azul suave

    # Linha radial do ponto inicial
    ax.plot([0, r0[0]], [0, r0[1]], [0, r0[2]], 'm--', linewidth=2, label='Linha Radial')  # Cor magenta

    # Marca o ponto inicial
    ax.scatter(y[0, 0], y[0, 1], y[0, 2], color='green', s=100, label='Ponto Inicial', edgecolors='black')

    # Marca o ponto final
    ax.scatter(y[-1, 0], y[-1, 1], y[-1, 2], color='blue', s=100, label='Ponto Final')
    ax.text(y[-1, 0], y[-1, 1], y[-1, 2], ' Fim', color='darkred', fontsize=12, weight='bold')

    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')
    ax.set_zlabel('Z (km)')
    ax.set_aspect('auto')

    # Define título e legenda
    ax.set_title('Visualização da Órbita', fontsize=16, fontweight='bold')
    ax.legend()

    # Ajustar os eixos para terem o mesmo tamanho
    set_axes_equal(ax)

    # Define uma perspectiva inicial para a visualização
    ax.view_init(elev=20., azim=30)

    plt.show()


def constroi_resultados(t, y, RE, r0, v0):
    y = np.asarray(y)
    if y.ndim == 1:
        raise ValueError("A matriz 'y' deve ser 2-dimensional.")

    r = np.linalg.norm(y[:, :3], axis=1)
    v = np.linalg.norm(y[:, 3:6], axis=1)
    rmax = np.max(r)
    rmin = np.min(r)
    imax = np.argmax(r)
    imin = np.argmin(r)
    v_at_rmax = v[imax]
    v_at_rmin = v[imin]
    altitude = r - RE

    imax_local, _ = find_peaks(altitude)
    imin_local, _ = find_peaks(-altitude)
    maxima = np.column_stack((t[imax_local], altitude[imax_local]))
    minima = np.column_stack((t[imin_local], altitude[imin_local]))

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
    plt.title('Altitude Over Time for Orbital Path', fontsize=14, fontweight='bold')
    plt.legend()
    plt.tight_layout()
    r_final = y[-1, :3]
    v_final = y[-1, 3:6]
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
    print("\nVetor posição no tempo final:", r_final)
    print("Vetor velocidade no tempo final:", v_final)

    plt.show()


