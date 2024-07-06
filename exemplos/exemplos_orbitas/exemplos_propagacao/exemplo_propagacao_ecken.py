import numpy as np

from src.domain.construtorderesultados.plots_orbitais import plota_orbita, plota_variacao_elementos_orbitais
from src.domain.modelos.orbitas.Orbita import Orbita
from src.domain.modelos.orbitas.propagacao.numerica.propagadores.propagacao_ecken import \
    propagacao_encke

r0 = np.array([-2.3845e3, 5.7290e3, 3.0505e3])
v0 = np.array([-7.3614, -2.9900, 1.6435, ])
mu = 398600.4418
raio_equatorial = 6378
orbita = Orbita.criar_pelo_vetor_de_estado(r0, v0, mu)
print(orbita.__repr__())
t, y = propagacao_encke(0, (24 * 60 * 60 * 2), orbita)
plota_orbita(y, raio_equatorial, r0)
plota_variacao_elementos_orbitais(t, y, orbita)
