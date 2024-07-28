import numpy as np

from src.construtorderesultados.plots_orbitais import plota_orbita
from src.domain.modelos.orbitas.Orbita import Orbita
from src.domain.modelos.orbitas.propagacao.numerica.propagadores.propagacao_numerica import propagacao_numerica
from src.domain.modelos.orbitas.utilidades.calculos_orbitais import calcular_periodo_orbital

ti = 0
r = np.array([8200, 8200, 0])
v = np.array([-4.5, 4.5, 0])

orbita = Orbita.criar_pelo_vetor_de_estado(r, v, 398600.4418)
tf = calcular_periodo_orbital(orbita.semi_eixo_maior, 398600.4418)

print("Periodo orbital:", tf)

t, y = propagacao_numerica(ti, tf, orbita)
plota_orbita(y, 6378, r)

