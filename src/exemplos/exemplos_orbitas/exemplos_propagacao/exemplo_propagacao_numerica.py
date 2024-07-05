from src.domain.utilidades_mecanica_orbital.Orbitas.coe_from_sv import coe_from_sv
from src.domain.utilidades_mecanica_orbital.Utilidades.calculos_orbitais import calcular_periodo_orbital
from src.domain.utilidades_mecanica_orbital.propagacao.numerica.propagadores.propaganumerica import propagacao_numerica
import numpy as np


r0 = [8200, 8200, 0]
v0 = [-4.5, 4.5, 1]
ti = 0
r = np.array([8200, 8200, 0])
v = np.array([-4.5, 4.5, 0])
a = coe_from_sv(r, v, 398600.4418)[0]
tf = calcular_periodo_orbital(a, 398600.4418)
print("Periodo orbital:", tf)
propagacao_numerica(ti, tf, r0, v0)
