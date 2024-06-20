import numpy as np

from src.domain.utilidades_mecanica_orbital.Orbitas.ModeloOrbita import Orbita

#
# Exemplo 5.2 do Tewari
r0 = np.array([-5000, 0, 12500])  # km
v0 = np.array([5, -8, 0])  # km/s
mu = 398600.4  # km^3/s^2 - Unidades coerentes com as entradas de posicao e velocidade

orb2 = Orbita.criar_pelo_vetor_de_estado(0, mu, r0, v0)
print('Exemplo 5.2 do Tewari')
print('Parametros orbitais identificados')
print(orb2.__repr__())
