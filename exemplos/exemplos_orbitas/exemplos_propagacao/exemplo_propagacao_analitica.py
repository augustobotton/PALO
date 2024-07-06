import matplotlib.pyplot as plt
import numpy as np

from src.domain.modelos.orbitas.ConstrutorOrbita import ConstrutorOrbita
from src.domain.modelos.orbitas.propagacao.analitica.propagacao import resolveEqKepler, matrizTransicaoEstado, \
    resolveEqKeplerHiperbolica
from src.domain.modelos.planeta.ConstrutorPlaneta import terra

# Constantes do problema
G = 6.67384e-11  # Constante de gravitação universal [m^3/kg*s^2]
ME = terra.massa  # Massa da Terra [kg]
RE = terra.raio_equatorial*1e3 # Raio equatorial da Terra [m]


# Condições iniciais
theta0 = 0  # Anomalia verdadeira inicial [rad] - no periastro
tau = 0  # Tempo de periastro [s] - supõe-se que ocorre em t=0

# Parâmetros orbitais para uma órbita elíptica
rp = 1.1 * RE  # Distância de periastro [m] - 10% maior que o raio equatorial
e = 0.3  # Excentricidade
a = rp / (1 - e)  # Semi-eixo maior [m]
p = (1 - e ** 2) * a  # Parâmetro orbital [m]

# Criação de uma órbita elíptica

orbita = ConstrutorOrbita().com_semi_eixo_maior(a).com_excentricidade(e).com_tempo_de_periastro(tau).construir()

# Cálculos orbitais
n = np.sqrt(terra.mut / a ** 3)
P = 2 * np.pi / n  # Período orbital [s]
N = 100  # Número de pontos
t = np.linspace(0, P, N)  # Vetor de tempo para um período completo

# Vetores posição e velocidade iniciais

h = np.sqrt(p * terra.mut)  # Quantidade de movimento angular específica [m^2/s]
r0 = p / (1 + e * np.cos(theta0))
R0 = np.array([r0 * np.cos(theta0), r0 * np.sin(theta0)])
V0 = np.array([-(terra.mut / h) * np.sin(theta0), (terra.mut / h) * (e + np.cos(theta0))])

# Inicialização de vetores para solução
E = np.zeros(N)  # Anomalia excêntrica
theta = np.zeros(N)  # Anomalia verdadeira
R = np.zeros((2, N))  # Posição no referencial perifocal
V = np.zeros((2, N))  # Velocidade no referencial perifocal

for i in range(N):
    E[i] = resolveEqKepler(t[i], orbita)
    theta[i] = 2 * np.arctan(np.sqrt((1 + e) / (1 - e)) * np.tan(E[i] / 2))
    theta[i] = theta[i] + 2 * np.pi if theta[i] < 0 else theta[i]

    PHI = matrizTransicaoEstado(theta[i], 0, orbita)
    R[:, i] = PHI[0, 0] * R0 + PHI[0, 1] * V0
    V[:, i] = PHI[1, 0] * R0 + PHI[1, 1] * V0

# Gráficos
plt.close("all")
plt.figure(figsize=(12, 8))

plt.subplot(2, 2, 1)
plt.plot(t, np.rad2deg(E), label="Anomalia Excêntrica")
plt.plot(t, np.rad2deg(theta), label="Anomalia Verdadeira")
plt.xlabel('Tempo [s]')
plt.ylabel('Ângulos [°]')
plt.grid(True)
plt.legend()

plt.subplot(2, 2, 3)
plt.plot(t, V[0, :] / 1e3, label="v_x")
plt.plot(t, V[1, :] / 1e3, label="v_y")
plt.xlabel('Tempo [s]')
plt.ylabel('Velocidades [km/s]')
plt.grid(True)
plt.legend()

plt.subplot(1, 2, 2)
xt = np.linspace(-RE, RE, N)
yt = np.sqrt(RE ** 2 - xt ** 2)
plt.plot(R[0, :] / 1e3, R[1, :] / 1e3, label="Órbita")
plt.plot(xt / 1e3, yt / 1e3, 'b', label="Terra")
plt.plot(xt / 1e3, -yt / 1e3, 'b')
plt.xlabel('X [km]')
plt.ylabel('Y [km]')
plt.axis('equal')
plt.grid(True)
plt.legend()

plt.tight_layout()
plt.show()




#####################################################################
# Constantes e parâmetros orbitais para uma órbita hiperbólica
e = 1.2  # Excentricidade representativa de uma hipérbole
rp = 1.1 * RE  # Distância de periastro [m] - 10% maior que o raio equatorial
a = -rp / (e - 1)  # Semi-eixo maior [m]
p = (1 - e ** 2) * a  # Parâmetro orbital [m]

# Criação da órbita hiperbólica

orbita2 = ConstrutorOrbita().com_semi_eixo_maior(a).com_excentricidade(e).com_tempo_de_periastro(0).construir()

# Determinação do 'período' baseado na distância de periastro para fins de simulação
n = np.sqrt(terra.mut / rp ** 3)
P = 2 * np.pi / n  # Período equivalente para cálculo de tempo [s]

# Vetor de tempo
N = 100  # Número de pontos para a simulação
t = np.linspace(0, P, N)

# Inicialização de vetores para posição e velocidade
r0 = p / (1 + e * np.cos(theta0))
R0 = np.array([r0 * np.cos(theta0), r0 * np.sin(theta0)])
V0 = np.array([-(terra.mut / h) * np.sin(theta0), (terra.mut / h) * (e + np.cos(theta0))])

# Vetores para armazenar a solução
H = np.zeros(N)  # Anomalia hiperbólica
theta = np.zeros(N)  # Anomalia verdadeira
R = np.zeros((2, N))  # Posição no referencial perifocal
V = np.zeros((2, N))  # Velocidade no referencial perifocal

for i in range(N):
    H[i] = resolveEqKeplerHiperbolica(t[i], orbita2)
    theta[i] = 2 * np.arctan(np.sqrt((e + 1) / (e - 1)) * np.tanh(H[i] / 2))
    theta[i] = theta[i] + 2 * np.pi if theta[i] < 0 else theta[i]

    PHI = matrizTransicaoEstado(theta[i], 0, orbita2)
    R[:, i] = PHI[0, 0] * R0 + PHI[0, 1] * V0
    V[:, i] = PHI[1, 0] * R0 + PHI[1, 1] * V0

# Gráficos
plt.figure(figsize=(12, 8))

plt.subplot(2, 2, 1)
plt.plot(t, np.rad2deg(H), label="Anomalia Hiperbólica (H)")
plt.plot(t, np.rad2deg(theta), label="Anomalia Verdadeira (θ)")
plt.xlabel('Tempo [s]')
plt.ylabel('Ângulos [°]')
plt.grid(True)
plt.legend()

plt.subplot(2, 2, 3)
plt.plot(t, V[0, :] / 1e3, label="v_x")
plt.plot(t, V[1, :] / 1e3, label="v_y")
plt.xlabel('Tempo [s]')
plt.ylabel('Velocidade [km/s]')
plt.grid(True)
plt.legend()

plt.subplot(1, 2, 2)
xt = np.linspace(-RE, RE, N)
yt = np.sqrt(RE ** 2 - xt ** 2)
plt.plot(R[0, :] / 1e3, R[1, :] / 1e3, label="Trajetória Hiperbólica")
plt.plot(xt / 1e3, yt / 1e3, 'b', label="Terra")
plt.plot(xt / 1e3, -yt / 1e3, 'b')
plt.xlabel('X [km]')
plt.ylabel('Y [km]')
plt.axis('equal')
plt.grid(True)
plt.legend()

plt.tight_layout()
plt.show()

