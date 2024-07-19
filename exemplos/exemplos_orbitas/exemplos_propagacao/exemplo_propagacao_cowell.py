import numpy as np

from src.domain.construtorderesultados.plots_orbitais import plota_orbita, constroi_resultados
from src.domain.modelos.orbitas.Orbita import Orbita
from src.domain.modelos.orbitas.manobras import delta_v_perigee_raise, aplicar_delta_v
from src.domain.modelos.orbitas.propagacao.numerica.propagadores.propagacao_cowell import propaga_cowell
from src.domain.modelos.orbitas.propagacao.numerica.propagadores.propagacao_numerica import propagacao_numerica
from src.domain.modelos.orbitas.utilidades.calculos_orbitais import calcular_periodo_orbital
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
criterio_de_parada = 200
args = [terra, 100, area, 2.2, criterio_de_parada]

# Tempo final da integração (120 dias em segundos)
tf = (24 * 60 * 60 * 120)

# Propagação da órbita utilizando o método de Cowell com perturbações de arrasto
t, y = propaga_cowell(0, tf, orbita, args)

# Plotagem da órbita
plota_orbita(y, terra.raio_equatorial, r0)

# Construção e exibição dos resultados
constroi_resultados(t, y, terra.raio_equatorial, r0, v0)

posicao_apos_120_dias = y[-1, :3]
velocidade_apos_120_dias = y[-1, 3:]

orbita_apos_120_dias = Orbita.criar_pelo_vetor_de_estado(posicao_apos_120_dias, velocidade_apos_120_dias, mu)
orbita_apos_120_dias.calcular_parametro_orbital()
deltav = delta_v_perigee_raise(orbita_apos_120_dias.calcula_periastro(),orbita_apos_120_dias.calcula_apoastro(),6577.9979+500 ,orbita_apos_120_dias.mu)

uma_orbita = calcular_periodo_orbital(orbita_apos_120_dias.semi_eixo_maior, orbita_apos_120_dias.mu)
t120, y120 = propagacao_numerica(0, uma_orbita*1.5, orbita_apos_120_dias)
plota_orbita(y120, terra.raio_equatorial, posicao_apos_120_dias)

print(f"Delta-V para aumentar a altitude perigeu em 500 km: {deltav} km/s")

pos = y120[-1, :3]
vel = y120[-1, 3:]


velocidadefinal = aplicar_delta_v(vel,deltav)
orbita_com_manobra = Orbita.criar_pelo_vetor_de_estado(pos, velocidadefinal, mu)
tfinal, yfinal = propagacao_numerica(0, uma_orbita, orbita_com_manobra)
constroi_resultados(tfinal, yfinal, terra.raio_equatorial, pos, velocidadefinal)
plota_orbita(yfinal, terra.raio_equatorial, pos)