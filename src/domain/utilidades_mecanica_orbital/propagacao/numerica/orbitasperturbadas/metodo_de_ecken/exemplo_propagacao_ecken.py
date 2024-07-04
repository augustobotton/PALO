import numpy as np
from src.domain.utilidades_mecanica_orbital.Orbitas.ConstrutorOrbita import ConstrutorOrbita
from src.domain.utilidades_mecanica_orbital.Orbitas.ModeloOrbita import Orbita
from src.domain.utilidades_mecanica_orbital.orbitalUtils.plotVariacoes import elementos_orbitais_resposta

from src.domain.utilidades_mecanica_orbital.propagacao.numerica.orbitasperturbadas.metodo_de_ecken.propagacao_ecken import \
    propagacao_encke
from matplotlib import pyplot as plt


def plot_orbit(y, R, r0):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    u, v = np.mgrid[0:2 * np.pi:100j, 0:np.pi:50j]
    x = R * np.cos(u) * np.sin(v)
    y_sphere = R * np.sin(u) * np.sin(v)
    z = R * np.cos(v)
    ax.plot_surface(x, y_sphere, z, color='grey', alpha=0.5)  # Desenha o planeta

    ax.plot(y[:, 0], y[:, 1], y[:, 2], 'k')  # Plota a trajetória da órbita
    ax.plot([0, r0[0]], [0, r0[1]], [0, r0[2]], 'r')  # Linha radial do ponto inicial
    ax.text(y[0, 0], y[0, 1], y[0, 2], 'o', color='red')  # Marca o ponto inicial
    ax.text(y[-1, 0], y[-1, 1], y[-1, 2], 'f', color='blue')  # Marca o ponto final

    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')
    ax.set_zlabel('Z (km)')
    ax.set_aspect('auto')
    plt.show()




r0 = np.array([ -2.3845e3,    5.7290e3,    3.0505e3])
v0 = np.array([-7.3614,   -2.9900,    1.6435,])
mu = 398600.4418
orbita = Orbita.criar_pelo_vetor_de_estado(r0, v0, mu)
print(orbita.__repr__())
t, y = propagacao_encke(0, (24*60*60*2), orbita)

plot_orbit(y, 6378, r0)

elementos_orbitais_resposta(y, t, orbita.mu)
