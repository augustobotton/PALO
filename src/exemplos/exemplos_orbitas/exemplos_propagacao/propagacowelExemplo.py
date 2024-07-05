import numpy as np

from src.domain.modelos.planeta.ConstrutorPlaneta import terra
from src.domain.utilidades_mecanica_orbital.Orbitas.ModeloOrbita import Orbita
from src.domain.utilidades_mecanica_orbital.Utilidades.plotVariacoes import elementos_orbitais_resposta
from src.domain.utilidades_mecanica_orbital.Utilidades.plotaOrbita import plot_orbit
from src.domain.utilidades_mecanica_orbital.propagacao.numerica.propagadores.propagaCowell import propaga_cowell


r0 = np.array([5872.94259799039,-662.622151725639,3007.48704280510])
v0 = np.array([-2.89355195466701,4.09603530867803,6.14446573555145])
mu = 398600.4418
orbita = Orbita.criar_pelo_vetor_de_estado(r0, v0, mu)
print(orbita.__repr__())
area = np.pi/4*(1**2)
args = [terra,100,area,2.2]
tf = (24 * 60 * 60 * 120)
t, y = propaga_cowell(0, tf, orbita, args)
plot_orbit(y, terra.raio_equatorial, r0)
