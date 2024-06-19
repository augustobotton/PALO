import numpy as np

from src.domain.utilidades_mecanica_orbital.Orbitas import exemplos_manobras
from src.domain.utilidades_mecanica_orbital.Orbitas.ModeloOrbita import Orbita

Re = 6378.137*1000  #[km]
mu = 3.986e14       #[m^2/s^2]
h_i = 250*1000      #[km]
r_i = Re + h_i
T_d = 86164         #[s]

r_f = np.cbrt(mu*(T_d/(2*np.pi))**2) #TODO add isso em calculos orbitais ou orbital utils

orbita_inicial = Orbita.circular(semi_eixo_maior=(h_i + Re),inclinacao=0)
d1,d2,t = exemplos_manobras.calcula_manobra_hohmann(orbita_inicial, r_f)

delta_v1 = np.sqrt(mu/r_i)*(np.sqrt((2*r_f)/(r_i+r_f))-1)
delta_v2 = np.sqrt(mu/r_f)*(1-np.sqrt((2*r_i)/(r_i+r_f)))
T_h = (np.sqrt((r_i+r_f)**3/(8*mu)))*np.pi

print(delta_v1)
print(d1)
print(delta_v2)
print(d2)
print(T_h )
print(t)
