import numpy as np

from src.domain.modelos.orbitas.Orbita import Orbita
from src.domain.modelos.orbitas.propagacao.numerica.propagadores.propagacao_cowell import propaga_cowell
from src.domain.modelos.orbitas.utilidades.plots import plot_orbit, constroe_resultados
from src.domain.modelos.planeta.ConstrutorPlaneta import terra

# Vetor de posição inicial em km
r0 = np.array([5872.94259799039, -662.622151725639, 3007.48704280510])
# Vetor de velocidade inicial em km/s
v0 = np.array([-2.89355195466701, 4.09603530867803, 6.14446573555145])
# Parâmetro gravitacional da Terra em km^3/s^2
mu = 398600.4418

# Criação de uma órbita a partir dos vetores de estado
orbita = Orbita.criar_pelo_vetor_de_estado(r0, v0, mu)
print(orbita.__repr__())

# Área de referência para o cálculo de arrasto (em m^2)
area = np.pi/4*(1**2)
# Argumentos para a dinâmica perturbada: planeta, massa (em kg), área de referência (em m^2), coeficiente de arrasto
args = [terra, 100, area, 2.2]

# Tempo final da integração (120 dias em segundos)
tf = (24 * 60 * 60 * 120)

# Propagação da órbita utilizando o método de Cowell com perturbações de arrasto
t, y = propaga_cowell(0, tf, orbita, args)

# Plotagem da órbita
plot_orbit(y, terra.raio_equatorial, r0)

# Construção e exibição dos resultados
constroe_resultados(terra.raio_equatorial, r0, t, v0, y)
