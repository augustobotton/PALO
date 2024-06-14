# Adaptado de André Luis da Silva
#
import numpy as np

from src.domain.utilidades_mecanica_orbital.orbitalUtils.determina_parametros_orbitais import \
	determina_parametros_orbitais

#
# Exemplo 5.2 do Tewari
r0 = np.array([-5000,0,12500]) # km
v0=np.array([5, -8,0]) # km/s
mu=398600.4 # km^3/s^2 - Unidades coerentes com as entradas de posicao e velocidade

orb=determina_parametros_orbitais(0, mu, r0, v0)
print('Exemplo 5.2 do Tewari')
print('Parametros orbitais identificados')
print('a = ', orb.semi_eixo_maior,'Km')
print('e = ',orb.excentricidade)
print('tau = ', orb.tempo_de_periastro,'s')
print('Omega = ', np.rad2deg(orb.raan),'°')
print('i = ', np.rad2deg(orb.inclinacao),'°')
print('omega = ',np.rad2deg(orb.arg_periastro),'°')