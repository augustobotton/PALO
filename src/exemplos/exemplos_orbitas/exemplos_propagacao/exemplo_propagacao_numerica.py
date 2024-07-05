import numpy as np

from src.domain.modelos.orbitas.propagacao.numerica.propagadores.propagacao_numerica import propagacao_numerica
from src.domain.modelos.orbitas.utilidades import calcular_periodo_orbital

r0 = [8200, 8200, 0]
v0 = [-4.5, 4.5, 1]
ti = 0
r = np.array([8200, 8200, 0])
v = np.array([-4.5, 4.5, 0])
a = coe_from_sv(r, v, 398600.4418)[0]
tf = calcular_periodo_orbital(a, 398600.4418)
print("Periodo orbital:", tf)
propagacao_numerica(ti, tf, r0, v0)
