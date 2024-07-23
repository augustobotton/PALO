import numpy as np

from src.domain.construtorderesultados.plota_ground_track import ground_track
from src.domain.construtorderesultados.plots_orbitais import plota_variacao_elementos_orbitais, plota_orbita
from src.domain.modelos.orbitas.Orbita import Orbita
from src.domain.modelos.orbitas.propagacao.numerica.propagadores.propagacao_encke import propagacao_encke
from src.domain.modelos.orbitas.propagacao.numerica.propagadores.propagacao_numerica import propagacao_numerica

posicao_inicial = np.array([-2384.46, 5729.01, 3050.46])
velocidade_inicial = np.array([-7.36138, -2.98997, 1.64354])
parametro_gravitacional = 398600

# Criação do objeto órbita inicial utilizando os vetores de estado
orbita = Orbita.criar_pelo_vetor_de_estado(posicao_inicial, velocidade_inicial, parametro_gravitacional)

# Impressão dos dados da órbita inicial
print(orbita.__repr__())

# Propagação da órbita utilizando o método de Encke
tempo, solucao_encke = propagacao_encke(0, (48*3600), orbita)
plota_orbita(solucao_encke, 6378.1370, posicao_inicial)
plota_variacao_elementos_orbitais(tempo, solucao_encke, orbita)
ground_track(tempo, solucao_encke)

# Propagação da órbita utilizando um método numérico
tempo2, solucao_numerica = propagacao_numerica(0, (48*3600), orbita)
plota_orbita(solucao_numerica, 6378.1370, posicao_inicial)
plota_variacao_elementos_orbitais(tempo2, solucao_numerica, orbita)
ground_track(tempo2, solucao_numerica)