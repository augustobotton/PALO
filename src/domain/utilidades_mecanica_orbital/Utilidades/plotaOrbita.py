import numpy as np
from matplotlib import pyplot as plt


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