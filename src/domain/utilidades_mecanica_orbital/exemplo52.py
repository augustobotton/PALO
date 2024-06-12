# Adaptado de André Luis da Silva
#
import numpy as np
from det_orbita import det_orbita
#
# Exemplo 5.2 do Tewari
r0 = np.array([-5000,0,12500]) # km
v0=np.array([5, -8,0]) # km/s
mu=398600.4 # km^3/s^2 - Unidades coerentes com as entradas de posicao e velocidade
par_orb=det_orbita(0,r0,v0,mu)
print('Exemplo 5.2 do Tewari')
print('Parametros orbitais identificados')
print('a = ', par_orb[0],'Km')
print('e = ',par_orb[1])
print('tau = ', par_orb[2],'s')
print('Omega = ', par_orb[3]*180/np.pi,'°')
print('i = ',par_orb[4]*180/np.pi,'°')
print('omega = ',par_orb[5]*180/np.pi,'°')