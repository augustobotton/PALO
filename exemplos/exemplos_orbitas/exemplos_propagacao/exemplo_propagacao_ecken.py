import numpy as np

from src.domain.construtorderesultados.plota_ground_track import ground_track
from src.domain.construtorderesultados.plots_orbitais import plota_variacao_elementos_orbitais, plota_orbita
from src.domain.modelos.orbitas.Orbita import Orbita
from src.domain.modelos.orbitas.propagacao.numerica.propagadores.propagacao_encke import propagacao_encke
from src.domain.modelos.orbitas.propagacao.numerica.propagadores.propagacao_numerica import propagacao_numerica

r0 = np.array([-2384.46, 5729.01, 3050.46])
v0 = np.array([-7.36138, -2.98997, 1.64354])
mu = 398600
orbita_inicial = Orbita.criar_pelo_vetor_de_estado(r0, v0, mu)
print(orbita_inicial.__repr__())
t, sol = propagacao_encke(0, (48*3600), orbita_inicial)


plota_orbita(sol,6378.1370, r0)
plota_variacao_elementos_orbitais(t, sol, orbita_inicial)
ground_track(t, sol)


r0 = np.array([-2384.46, 5729.01, 3050.46])
v0 = np.array([-7.36138, -2.98997, 1.64354])
mu = 398600
orbita_inicial = Orbita.criar_pelo_vetor_de_estado(r0, v0, mu)
print(orbita_inicial.__repr__())
t, sol = propagacao_numerica(0, (48*3600), orbita_inicial)

plota_orbita(sol,6378.1370, r0)
plota_variacao_elementos_orbitais(t, sol, orbita_inicial)
ground_track(t, sol)



t
