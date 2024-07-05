import numpy as np
from src.domain.utilidades_mecanica_orbital.Orbitas.ModeloOrbita import Orbita
from src.domain.utilidades_mecanica_orbital.Utilidades.plotVariacoes import elementos_orbitais_resposta
from src.domain.utilidades_mecanica_orbital.Utilidades.plotaOrbita import plot_orbit
from src.domain.utilidades_mecanica_orbital.propagacao.numerica.propagadores.metodo_de_ecken.propagacao_ecken import \
    propagacao_encke

r0 = np.array([-2.3845e3, 5.7290e3, 3.0505e3])
v0 = np.array([-7.3614, -2.9900, 1.6435, ])
mu = 398600.4418
orbita = Orbita.criar_pelo_vetor_de_estado(r0, v0, mu)
print(orbita.__repr__())
t, y = propagacao_encke(0, (24 * 60 * 60 * 2), orbita)
plot_orbit(y, 6378, r0)
elementos_orbitais_resposta(y, t, orbita.mu)
