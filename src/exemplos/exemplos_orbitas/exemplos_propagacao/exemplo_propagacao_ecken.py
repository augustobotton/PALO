import numpy as np

from src.domain.modelos.orbitas import Orbita
from src.domain.modelos.orbitas.propagacao.numerica.propagadores.propagacao_ecken import \
    propagacao_encke
from src.domain.modelos.orbitas.utilidades import elementos_orbitais_resposta
from src.domain.modelos.orbitas.utilidades import plot_orbit

r0 = np.array([-2.3845e3, 5.7290e3, 3.0505e3])
v0 = np.array([-7.3614, -2.9900, 1.6435, ])
mu = 398600.4418
orbita = Orbita.criar_pelo_vetor_de_estado(r0, v0, mu)
print(orbita.__repr__())
t, y = propagacao_encke(0, (24 * 60 * 60 * 2), orbita)
plot_orbit(y, 6378, r0)
elementos_orbitais_resposta(y, t, orbita.mu)
