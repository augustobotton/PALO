# Exemplo de uso com 3 estágios
import numpy as np

from src.domain.modelos.orbitas.Orbita import Orbita
from src.domain.modelos.orbitas.utilidades import calculos_orbitais
from src.domain.modelos.planeta.ConstrutorPlaneta import terra

#Inicialmente é calculado a velocidade orbital necessária para a órbita desejada
orbita_alvo = Orbita.circular(1000.164e3 + terra.raio_equatorial, np.deg2rad(40))

velocidade_orbital_desejada = calculos_orbitais.calcula_velocidade_orbital(terra.mut, terra.raio_equatorial+orbita_alvo.semi_eixo_maior,orbita_alvo.semi_eixo_maior)
print(velocidade_orbital_desejada)
impulso_velocidade = 9.5 * 1000  # Impulso de velocidade requerido para órbita baixa (em m/s)
razoes_estruturais = [0.07, 0.05, 0.05]  # Razões estruturais
velocidade_exaustao1 = 200 * 9.81  # Velocidade de exaustão do primeiro estágio (em m/s)
exaustao_normalizada = [1, 1.5, 1.75]  # Velocidades de exaustão normalizadas
velocidades_exaustao = [velocidade_exaustao1 * b for b in exaustao_normalizada]
carga_util = 5000  # Carga útil (em kg)
numero_estagios = 3
estagios_paralelo = None  # Sem estágios em paralelo







