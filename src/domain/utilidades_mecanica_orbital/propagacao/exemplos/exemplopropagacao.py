import matplotlib.pyplot as plt
import numpy as np

from src.domain.utilidades_mecanica_orbital.propagacao.propagaOrbKepleriana import *

#
# Script com exemplos de propagacao de orbita
#
# Constantes gerais do problema
G = 6.67384e-11  # [m^3/kg*s^2] Constante de gravitacao universal
Me = 5.9722e24  # [kg] Massa da Terra
Re = 6.37814e6  # [m] Raio equatorial da Terra
mu = G * Me  # Constante gravitacional da Terra
# Condicoes iniciais dos exemplos
theta0 = 0  # [rad] Anomalia verdadeira inicial - supor o periastro
tau = 0  # [s] Tempo de  periastro - supor que ocorre em t=0
# Exportacao para outras funcoes usando um arquivo de parametros

#


# Aplicacao em orbita eliptica
#
rp = 1.1 * Re  # [m] - Distancia de periastro - valor generico 10% maior que Re
e = 0.3  # Excentricidade representativa de uma elipse
a = rp / (1 - e)  # [m] Semi eixo maior
p = (1 - e ** 2) * a  # [m] Parametro

criar_orbita = ConstrutorOrbita()
orbita = criar_orbita.com_semi_eixo_maior(a).com_excentricidade(e).com_tempo_de_periastro(0).construir()
#############################################################################
# Cálculos - solução do exemplo de órbita eliptica
# Período
n = np.sqrt(mu / a ** 3);
f = n / (2 * np.pi);
P = 1 / f
# Vetor de tempo
N = 100;
t = np.linspace(0, P, N)  # Vetor para um período de orbita, 100 divisões
# Vetores posição e velocidade inicial, no referencial perifocal, escritos
# em coordenadas retangulares
h = np.sqrt(p * mu)  # [m^2/s] Quantidade de movimento angular especifica
r0 = p / (1 + e * np.cos(theta0))
R0 = np.array([r0 * np.cos(theta0), r0 * np.sin(theta0)])
V0 = np.array([-(mu / h) * np.sin(theta0), (mu / h) * (e + np.cos(theta0))])
# Matrizes que guardarão a solução
E = np.zeros(N)  # Anomalia excentrica
theta = np.zeros(N)  # Anomalia verdadeira
R = np.zeros((2, N))  # Vetor posicao no referencial perifocal em coordenadas retangulares
V = np.zeros((2, N))  # Vetor velocidade no referencial perifocal em coordenadas retangulares
# Cálculo da órbita em cada instante de tempo
for i in range(N):
    # Resolve a equacao de Kepler, determinando a anomalia excêntrica em
    # cada instante de tempo
    E[i] = resolveEqKepler(t[i], orbita)
    # Determina a anomalia verdadeira a partir da excêntrica
    theta[i] = 2 * np.arctan(np.sqrt((1 + e) / (1 - e)) * np.tan(E[i] / 2))
    if theta[i] < 0:
        theta[i] = theta[i] + 2 * np.pi
    # Determina a matriz de transição de estado a partir de theta
    PHI = matrizTransicaoEstado(theta[i], 0, orbita)
    # Determina posição e velocidade perifocal em função da anomalia
    # verdadeira pela matriz de transição de estado
    R[:, i] = PHI[0, 0] * R0 + PHI[0, 1] * V0
    V[:, i] = PHI[1, 0] * R0 + PHI[1, 1] * V0

# Ilustração da Terra
xt = np.linspace(-Re, Re, N)
yt = np.sqrt(Re ** 2 - xt ** 2)
# Graficos
plt.close("all")
#
plt.figure(1);
#
plt.subplot(2, 2, 1)
plt.plot(t, E * 180 / np.pi, label="E")
plt.plot(t, theta * 180 / np.pi, label="theta")
plt.xlabel('t [s]')
plt.ylabel('theta, E [°]')
plt.grid()
plt.legend()
#
plt.subplot(2, 2, 3)
plt.plot(t, V[0, :] / 1e3, label="v_x")
plt.plot(t, V[1, :] / 1e3, label="v_y")
plt.xlabel('t [s]');
plt.ylabel('v_x, v_y [km/s]')
plt.grid();
plt.legend()
#
plt.subplot(1, 2, 2);
plt.plot(R[0, :] / 1e3, R[1, :] / 1e3)
plt.plot(xt / 1e3, yt / 1e3, 'b')
plt.plot(xt / 1e3, -yt / 1e3, 'b')
plt.xlabel('X [km]')
plt.ylabel('Y [km]')
plt.axis('equal')
plt.grid()
plt.legend()

##################################################
# Exemplo de órbita hiperbólica
rp = 1.1 * Re  # [m] Distancia de periastro - valor generico 10% maior que Re
e = 1.2  # Excentricidade representativa de uma hipérbole
a = -rp / (e - 1)  # Semi eixo maior
p = (1 - e ** 2) * a  # Parametro
orbita2 = criar_orbita.com_semi_eixo_maior(a).com_excentricidade(e).com_tempo_de_periastro(0).construir()
# Periodo de um orbita circular de mesma distancia de periastro
n = np.sqrt(mu / rp ** 3)
f = n / (2 * np.pi)
P = 1 / f

# Vetor de tempo
N = 100
t = np.linspace(0, P, N)  # Vetor para um período de orbita, 100 divisões
# Vetores posicao e velocidade inicial, no referencial perifocal, escritos
# em coordenadas retangulares
r0 = p / (1 + e * np.cos(theta0))
R0 = np.array([r0 * np.cos(theta0), r0 * np.sin(theta0)])
V0 = np.array([-(mu / h) * np.sin(theta0), (mu / h) * (e + np.cos(theta0))])
# Matrizes que guardarao a solucao
H = np.zeros(N)  # Anomalia hiperbólica
theta = np.zeros(N)  # Anomalia verdadeira
R = np.zeros((2, N))  # Vetor posicao no referencial perifocal em coordenadas retangulares
V = np.zeros((2, N))  # Vetor velocidade no referencial perifocal em coordenadas retangulares

# Cálculo da órbita em cada instante de tempo
for i in range(N):
    # Resolve a equacao de Kepler hiperbolica, determinando a anomalia hiperbolica em
    # cada tempo
    H[i] = resolveEqKeplerHiperbolica(t[i], orbita2)
    # Determina a anomalia verdadeira a partir da hiperbolica
    theta[i] = 2 * np.arctan(np.sqrt((e + 1) / (e - 1)) * np.tanh(H[i] / 2))
    # Determina a matriz de transicao de estado a partir de theta
    PHI = matrizTransicaoEstado(theta[i], 0,orbita)
    # Determina posicao e velocidade perifocal em funcao de suas condicoes
    # iniciais pela matriz de transicao de estado
    R[:, i] = PHI[0, 0] * R0 + PHI[0, 1] * V0
    V[:, i] = PHI[1, 0] * R0 + PHI[1, 1] * V0

# Graficos
plt.figure(2);
#
plt.subplot(2, 2, 1)
plt.plot(t, H * 180 / np.pi, label="H")
plt.plot(t, theta * 180 / np.pi, label="theta")
plt.xlabel('t [s]')
plt.ylabel('theta, H [°]')
plt.grid()
plt.legend()
#
plt.subplot(2, 2, 3)
plt.plot(t, V[0, :] / 1e3, label="v_x")
plt.plot(t, V[1, :] / 1e3, label="v_y")
plt.xlabel('t [s]')
plt.ylabel('v_x, v_y [km/s]')
plt.grid()
plt.legend()
#
plt.subplot(1, 2, 2)
plt.plot(R[0, :] / 1e3, R[1, :] / 1e3)
plt.plot(xt / 1e3, yt / 1e3, 'b')
plt.plot(xt / 1e3, -yt / 1e3, 'b')
plt.xlabel('X [km]')
plt.ylabel('Y [km]')
plt.axis('equal')
plt.grid()
plt.legend()

#################################################

# Aplicacao em orbita parabolica
rp = 1.1 * Re  # [m] Distancia de periastro - valor generico 10% maior que Re
p = rp * 2  # Parametro
e = 1


# Exportacao para outras funcoes usando um arquivo de parametros
orbita3 = criar_orbita.com_parametro(p).com_excentricidade(e).com_tempo_de_periastro(0).construir()
# Periodo de um orbita circular de mesma distancia de periastro
n = np.sqrt(mu / rp ** 3)
f = n / (2 * np.pi)
P = 1 / f
# Vetor de tempo
N = 100
t = np.linspace(0, P, N)  # Vetor para um período de orbita, 100 divisões
# Vetores posicao e velocidade inicial, no referencial perifocal, escritos
# em coordenadas retangulares
r0 = p / (1 + e * np.cos(theta0))
R0 = np.array([r0 * np.cos(theta0), r0 * np.sin(theta0)])
V0 = np.array([-(mu / h) * np.sin(theta0), (mu / h) * (e + np.cos(theta0))])
# Matrizes que guardarao a solucao
theta = np.zeros(N)  # Anomalia verdadeira
R = np.zeros((2, N))  # Vetor posicao no referencial perifocal em coordenadas retangulares
V = np.zeros((2, N))  # Vetor velocidade no referencial perifocal em coordenadas retangulares
# Cálculo da órbita em cada instante de tempo
for i in range(N):
    # Resolve a equacao de Kepler de Barker, determinando a anomalia verdadeira em
    # cada tempo
    theta[i] = resolveEqBarker(t[i], orbita3)
    # Determina a matriz de transicao de estado a partir de theta
    PHI = matrizTransicaoEstado(theta[i],0, orbita3)
    # Determina posicao e velocidade perifocal em funcao de suas condicoes
    # iniciais pela matriz de transicao de estado
    R[:, i] = PHI[0, 0] * R0 + PHI[0, 1] * V0
    V[:, i] = PHI[1, 0] * R0 + PHI[1, 1] * V0

# Gráficos
plt.figure(3)
#
plt.subplot(2, 2, 1)
plt.plot(t, theta * 180 / np.pi)
plt.xlabel('t [s]')
plt.ylabel('theta [°]')
plt.grid()
plt.legend()
#
plt.subplot(2, 2, 3)
plt.plot(t, V[0, :] / 1e3, label="v_x")
plt.plot(t, V[1, :] / 1e3, label="v_y")
plt.xlabel('t [s]')
plt.ylabel('v_x, v_y [km/s]')
plt.grid()
plt.legend()
#
plt.subplot(1, 2, 2)
plt.plot(R[0, :] / 1e3, R[1, :] / 1e3)
plt.plot(xt / 1e3, yt / 1e3, 'b');
plt.plot(xt / 1e3, -yt / 1e3, 'b')
plt.xlabel('X [km]')
plt.ylabel('Y [km]')
plt.axis('equal')
plt.grid()
plt.legend()
plt.show()
