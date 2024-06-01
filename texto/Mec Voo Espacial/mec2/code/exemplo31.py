import parametros
import numpy as np
import matplotlib.pyplot as plt
from gravidade_axissim import grav_axissimetrico
# Exemplo 3.1 do Tewari. Comparar 2 modelos de gravidade Terrestre ao longo da seguinte
# trajetoria: delta=h-100 (0<=h<=200 km), delta eh a latitude em graus.
# Modelos comparados: Terra esferica e Terra axissimetrica.
#
# Definicao das constantes do problema e insercao no modulo de variaveis globais
Re=6378.14e3 # m - Raio equatorial do planeta
parametros.RE=Re
G=6.67259e-11 # m^3/kg/s^2 - Constante de gravitacao universal no SI
parametros.G=G
M=5.972e24 # kg - Massa do planeta
parametros.M=M
J2=0.00108263;J3=-0.00000254;J4=-0.00000161 # Constantes de Jeffery
parametros.J2=J2;parametros.J3=J3;parametros.J4=J4
#
# Resolucao do exercicio
#
# Define um vetor de altitude
N=100
h=np.linspace(0,200,N) # km
# Latitude na trajetoria
delta=h-100 # graus
# O modelo gravitacional de corpo axissimetrico eh funcao da distancia radial
# e da colatitude. Abaixo sao feitas as conversoes necessarias.
# Calculo da distancia radial. Assume-se que a altitude eh com respeito ao raio equatorial
r=parametros.RE+h*1000 # m - A altitude deve ser convertida para metros

# Colatitude
phi=np.pi/2-delta*np.pi/180 # rad - a latitude deve ser convertida para radianos
#


# Calculo da gravidade comparando os dois modelos
g=np.empty(N) # Inicializa vetor para salvar a gravidade de modelo esferico
gr=np.empty(N) # Inicializa o vetor para salvar a componente radial do modelo axissimetrico
gphi=np.empty(N) # Inicializa o vetor para salvar a componente colatitudinal (sul) do modelo axissime

for i in range(N):
    # Modelo de planeta esferico
    g[i]=-G*M/r[i]**2;
    # Modelo de planeta axissimetrico
    gr[i],gphi[i]=grav_axissimetrico(r[i],phi[i])

# Graficos
plt.close('all')
plt.figure(1)
plt.subplot(2,1,1);plt.plot(delta,g, label = "esferico")
plt.plot(delta,gr, label = "axissimetrico")
plt.xlabel('latitude - graus');plt.ylabel('gravidade m/s^2')
plt.grid();plt.legend();
plt.subplot(212);plt.plot(delta,gphi, label = "axissimetrico")
plt.xlabel('latitude - graus');plt.ylabel('gravidade m/s^2')
plt.grid();plt.legend();plt.show()