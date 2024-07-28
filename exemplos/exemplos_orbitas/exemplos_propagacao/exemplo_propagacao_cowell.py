import numpy as np

from src.construtorderesultados.plots_orbitais import plota_orbita
from src.domain.modelos.orbitas.Orbita import Orbita
from src.domain.modelos.orbitas.manobras import calcula_delta_v_eleva_perigeu, aplicar_delta_v
from src.domain.modelos.orbitas.propagacao.numerica.propagadores.propagacao_cowell import propaga_cowell
from src.domain.modelos.orbitas.propagacao.numerica.propagadores.propagacao_numerica import propagacao_numerica
from src.domain.modelos.orbitas.utilidades.calculos_orbitais import calcular_periodo_orbital
from src.domain.modelos.orbitas.utilidades.outras_utilidades import constroi_resultados
from src.domain.modelos.planeta.ConstrutorPlaneta import terra

posicao_inicial = np.array([5872.94259799039, -662.622151725639, 3007.48704280510])
velocidade_inicial = np.array([-2.89355195466701, 4.09603530867803, 6.14446573555145])

# Parâmetro gravitacional da Terra em km^3/s^2
parametro_gravitacional = 398600.4418

# Criação de uma órbita a partir dos vetores de estado
orbita_inicial = Orbita.criar_pelo_vetor_de_estado(posicao_inicial, velocidade_inicial, parametro_gravitacional)
print(orbita_inicial.__repr__())

# Área de referência para o cálculo de arrasto (em m^2)
area_referencia = np.pi/4*(1**2)

# Argumentos para a dinâmica perturbada: planeta, massa (em kg), área de referência (em m^2), coeficiente de arrasto, critério de parada

altitude_critica = 200 # Indica para o solver para parar a integração
CD = 2.2 # Coeficiente de arrasto
massa = 100 #  kg

parametros_dinamica = [terra, massa, area_referencia, CD, altitude_critica]

# Tempo final da integração (120 dias em segundos)
tempo_final = (24 * 60 * 60 * 120)

# Propagação da órbita utilizando o método de Cowell com perturbações de arrasto
tempo_propagacao, estado_propagado = propaga_cowell(0, tempo_final, orbita_inicial, parametros_dinamica)

# Plotagem da órbita propagada considerando perturbações
plota_orbita(estado_propagado, terra.raio_equatorial, posicao_inicial)

# Construção e exibição dos resultados da propagação
constroi_resultados(tempo_propagacao, estado_propagado, terra.raio_equatorial, posicao_inicial, velocidade_inicial)

# Obtendo a posição e velocidade após 120 dias de propagação
posicao_apos_120_dias = estado_propagado[-1, :3]
velocidade_apos_120_dias = estado_propagado[-1, 3:]

# Criação de uma nova órbita a partir da posição e velocidade após 120 dias
orbita_apos_120_dias = Orbita.criar_pelo_vetor_de_estado(posicao_apos_120_dias, velocidade_apos_120_dias, parametro_gravitacional)

# Cálculo do Delta-V necessário para aumentar a altitude do perigeu em 500 km
delta_v = calcula_delta_v_eleva_perigeu(orbita_apos_120_dias.calcula_periastro(), orbita_apos_120_dias.calcula_apoastro(), 6577.9979 + 150 , orbita_apos_120_dias.mu)

# Cálculo do período orbital da nova órbita, usado para simular apenas uma órbita em torno da terra
periodo_orbital = calcular_periodo_orbital(orbita_apos_120_dias.semi_eixo_maior, orbita_apos_120_dias.mu)

# Propagação da nova órbita por um período e meio
tempo_propagacao_120_dias, estado_apos_120_dias = propagacao_numerica(0, periodo_orbital * 1.5, orbita_apos_120_dias)
plota_orbita(estado_apos_120_dias, terra.raio_equatorial, posicao_apos_120_dias)



# Obtendo a posição e velocidade ao final da propagação
posicao_final = estado_apos_120_dias[-1, :3]
velocidade_final = estado_apos_120_dias[-1, 3:]

# Exibição do Delta-V calculado
print(f"Delta-V para aumentar a altitude perigeu em 150 km: {delta_v} km/s")

nova_velocidade = aplicar_delta_v(velocidade_final, delta_v)

# Criação de uma nova órbita considerando a manobra Delta-V
orbita_corrigida = Orbita.criar_pelo_vetor_de_estado(posicao_final, nova_velocidade, parametro_gravitacional)

# Propagação da órbita com a manobra Delta-V utilizando o método de Cowell
tempo_final_corrigido, estado_corrigido = propaga_cowell(0, tempo_final, orbita_corrigida, parametros_dinamica)
constroi_resultados(tempo_final_corrigido, estado_corrigido, terra.raio_equatorial, posicao_final, nova_velocidade)
plota_orbita(estado_corrigido, terra.raio_equatorial, posicao_final)
