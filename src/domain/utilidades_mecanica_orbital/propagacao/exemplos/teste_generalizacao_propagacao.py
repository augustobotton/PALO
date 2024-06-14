import matplotlib.pyplot as plt
import numpy as np
from src.domain.utilidades_mecanica_orbital.propagacao.propagaOrbKepleriana import *

# Constantes gerais do problema
G = 6.67384e-11  # [m^3/kg*s^2] Constante de gravitação universal
Me = 5.9722e24  # [kg] Massa da Terra
Re = 6.37814e6  # [m] Raio equatorial da Terra
mu = G * Me  # Constante gravitacional da Terra


# Função para calcular parâmetros orbitais
def calcular_parametros_orbitais(a, e):
    if a > 0:
        p = (1 - e ** 2) * a  # [m] Parâmetro
        n = np.sqrt(mu / a ** 3)
        f = n / (2 * np.pi)
        P = 1 / f
        return p, P
    else:
        return None, None


# Função para resolver órbita elíptica
def resolver_orbita_eliptica(a, e, theta0, tau):
    p, P = calcular_parametros_orbitais(a, e)
    if p is None or P is None:
        print("Erro: Parâmetros orbitais inválidos.")
        return

    orbita = ConstrutorOrbita().com_semi_eixo_maior(a).com_excentricidade(e).com_tempo_de_periastro(tau).construir()
    t = np.linspace(0, P, 100)
    h = np.sqrt(p * mu)
    r0 = p / (1 + e * np.cos(theta0))
    R0 = np.array([r0 * np.cos(theta0), r0 * np.sin(theta0)])
    V0 = np.array([-(mu / h) * np.sin(theta0), (mu / h) * (e + np.cos(theta0))])
    E, theta, R, V = calcular_orbita(t, orbita, e, R0, V0)
    plotar_orbita(t, E, theta, R, V, Re, 1, tipo='eliptica')


# Função para resolver órbita hiperbólica
def resolver_orbita_hiperbolica(a, e, theta0, tau):
    p, P = calcular_parametros_orbitais(a, e)
    if p is None or P is None:
        print("Erro: Parâmetros orbitais inválidos.")
        return

    orbita = ConstrutorOrbita().com_semi_eixo_maior(a).com_excentricidade(e).com_tempo_de_periastro(tau).construir()
    t = np.linspace(0, P, 100)
    h = np.sqrt(p * mu)
    r0 = p / (1 + e * np.cos(theta0))
    R0 = np.array([r0 * np.cos(theta0), r0 * np.sin(theta0)])
    V0 = np.array([-(mu / h) * np.sin(theta0), (mu / h) * (e + np.cos(theta0))])
    H, theta, R, V = calcular_orbita(t, orbita, e, R0, V0, tipo='hiperbolica')
    plotar_orbita(t, H, theta, R, V, Re, 2, tipo='hiperbolica')


# Função para resolver órbita parabólica
def resolver_orbita_parabolica(p, e, theta0, tau):
    orbita = ConstrutorOrbita().com_parametro(p).com_excentricidade(e).com_tempo_de_periastro(tau).construir()
    n = np.sqrt(mu / p ** 3)
    f = n / (2 * np.pi)
    P = 1 / f
    t = np.linspace(0, P, 100)
    h = np.sqrt(p * mu)
    r0 = p / (1 + e * np.cos(theta0))
    R0 = np.array([r0 * np.cos(theta0), r0 * np.sin(theta0)])
    V0 = np.array([-(mu / h) * np.sin(theta0), (mu / h) * (e + np.cos(theta0))])
    theta, R, V = calcular_orbita_parabolica(t, orbita, R0, V0)
    plotar_orbita(t, theta, R, V, Re, 3, tipo='parabolica')


# Função para calcular a órbita
def calcular_orbita(t, orbita, e, R0, V0, tipo='eliptica'):
    N = len(t)
    E, theta = np.zeros(N), np.zeros(N)
    R, V = np.zeros((2, N)), np.zeros((2, N))
    H = np.zeros(N)  # Adicionando a definição de H
    for i in range(N):
        if tipo == 'eliptica':
            E[i] = resolveEqKepler(t[i], orbita)
            theta[i] = 2 * np.arctan(np.sqrt((1 + e) / (1 - e)) * np.tan(E[i] / 2))
        else:
            H[i] = resolveEqKeplerHiperbolica(t[i], orbita)
            theta[i] = 2 * np.arctan(np.sqrt((e + 1) / (e - 1)) * np.tanh(H[i] / 2))
        if theta[i] < 0:
            theta[i] += 2 * np.pi
        PHI = matrizTransicaoEstado(theta[i], 0, orbita)
        R[:, i] = PHI[0, 0] * R0 + PHI[0, 1] * V0
        V[:, i] = PHI[1, 0] * R0 + PHI[1, 1] * V0
    return E if tipo == 'eliptica' else H, theta, R, V


# Função para calcular órbita parabólica
def calcular_orbita_parabolica(t, orbita, R0, V0):
    N = len(t)
    theta, R, V = np.zeros(N), np.zeros((2, N)), np.zeros((2, N))
    for i in range(N):
        theta[i] = resolveEqBarker(t[i], orbita)
        PHI = matrizTransicaoEstado(theta[i], 0, orbita)
        R[:, i] = PHI[0, 0] * R0 + PHI[0, 1] * V0
        V[:, i] = PHI[1, 0] * R0 + PHI[1, 1] * V0
    return theta, R, V


# Função para plotar as órbitas
def plotar_orbita(t, E, theta, R, V, Re, figura, tipo='eliptica'):
    plt.figure(figura)
    plt.subplot(2, 2, 1)
    if tipo == 'eliptica':
        plt.plot(t, E * 180 / np.pi, label="E")
    else:
        plt.plot(t, E * 180 / np.pi, label="H")
    plt.plot(t, theta * 180 / np.pi, label="theta")
    plt.xlabel('t [s]')
    plt.ylabel('theta, E/H [°]')
    plt.grid()
    plt.legend()

    plt.subplot(2, 2, 3)
    plt.plot(t, V[0, :] / 1e3, label="v_x")
    plt.plot(t, V[1, :] / 1e3, label="v_y")
    plt.xlabel('t [s]')
    plt.ylabel('v_x, v_y [km/s]')
    plt.grid()
    plt.legend()

    xt = np.linspace(-Re, Re, len(t))
    yt = np.sqrt(Re ** 2 - xt ** 2)

    plt.subplot(1, 2, 2)
    plt.plot(R[0, :] / 1e3, R[1, :] / 1e3, label="Orbita")
    plt.plot(xt / 1e3, yt / 1e3, 'b', label="Terra")
    plt.plot(xt / 1e3, -yt / 1e3, 'b')
    plt.xlabel('X [km]')
    plt.ylabel('Y [km]')
    plt.axis('equal')
    plt.grid()
    plt.legend()
    plt.show()


# Parâmetros iniciais
theta0 = 0  # [rad] Anomalia verdadeira inicial
tau = 0  # [s] Tempo de periastro

# Resolver órbitas
resolver_orbita_eliptica(1.1 * Re / (1 - 0.3), 0.3, theta0, tau)
resolver_orbita_hiperbolica(-1.1 * Re / (1.2 - 1), 1.2, theta0, tau)
resolver_orbita_parabolica(2 * 1.1 * Re, 1, theta0, tau)
