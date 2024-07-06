import matplotlib.pyplot as plt
import numpy as np

from src.domain.modelos.planeta.ConstrutorPlaneta import terra
from src.domain.modelos.planeta.ModeloGravidadePlanetaAxisimetrico import calcular_gravidade_axisimetrico

# Constante de gravitação universal (m^3/kg/s^2)
G = 6.67259e-11

# Ajustando as unidades de constantes do planeta
terra.mut *= 1e3 # Converte unidade devido a incopatibilidade dos modelos de órbita e lançamento
terra.raio_equatorial *= 1e3

# Definindo o vetor de altitudes
num_pontos = 100
altitudes_km = np.linspace(0, 200, num_pontos)  # km

# Convertendo latitude para radianos
latitudes_rad = np.deg2rad(altitudes_km  - 100)

# Calculando a distância radial (em metros)
distancias_radiais = terra.raio_equatorial + (altitudes_km * 1e3)

# Preparando vetores para armazenar resultados gravitacionais
gravidade_esferica = np.empty(num_pontos)
gravidade_radial = np.empty(num_pontos)
gravidade_colatitudinal = np.empty(num_pontos)

for i in range(num_pontos):
    # Calculando para modelo de planeta esférico
    gravidade_esferica[i] = (-G * terra.massa) / distancias_radiais[i] ** 2

    # Calculando para modelo de planeta axissimétrico
    gravidade_radial[i], gravidade_colatitudinal[i] = calcular_gravidade_axisimetrico(distancias_radiais[i],
                                                                                      latitudes_rad[i], terra)

# Convertendo radianos para graus para uso no gráfico
latitudes_graus = np.rad2deg(latitudes_rad)

# Criando gráficos comparativos
plt.figure(figsize=(10, 8))

# Gráfico para a gravidade esférica e axissimétrica radial
plt.subplot(2, 1, 1)
plt.plot(latitudes_graus, gravidade_esferica, label="Gravidade Esférica")
plt.plot(latitudes_graus, -gravidade_radial*1e6, label="Gravidade Axissimétrica Radial")
plt.xlabel('Latitude (graus)')
plt.ylabel('Gravidade (m/s²)')
plt.title('Comparação da Gravidade Esférica e Axissimétrica Radial')
plt.grid(True)
plt.legend()

# Gráfico para a gravidade axissimétrica colatitudinal
plt.subplot(2, 1, 2)
plt.plot(latitudes_graus, -gravidade_colatitudinal*1e6, label="Gravidade Axissimétrica Colatitudinal")
plt.xlabel('Latitude (graus)')
plt.ylabel('Gravidade (m/s²)')
plt.title('Gravidade Axissimétrica Colatitudinal')
plt.grid(True)
plt.legend()

plt.tight_layout()
plt.show()
