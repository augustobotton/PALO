import pickle

import numpy as np

from src.domain.modelos.Simulacao import Simulacao
from src.domain.modelos.foguete.aerodinamica.ModeloAerodinamico import ModeloAerodinamico
from src.domain.modelos.foguete.estrutura.ConstrutorModeloEstrutural import ConstrutorModeloEstrutural
from src.domain.modelos.foguete.propulsao.ConstrutorModeloPropulsivo import ConstrutorModeloPropulsivo
from src.domain.modelos.foguete.veiculo_lancador.ConstrutorFoguete import ConstrutorDeFoguete
from src.domain.modelos.orbitas.Orbita import Orbita
from src.domain.modelos.planeta.ConstrutorBaseDeLancamento import ConstrutorBaseDeLancamento
from src.domain.modelos.planeta.ConstrutorPlaneta import ConstrutorPlaneta
from src.domain.modelos.planeta.ModeloAtmosferico import ModeloAtmosferico

# Estágio 1
Isp_1 = 290  # s
F_1 = 33400  # kN
m0_1 = 2780000  # kg
mp_1 = 1997000  # kg
me_1 = 106000  # kg
epsilon_1 = 0.050
lambda_1 = 0.321
burn_time_1 = 168  # s
area_cross_1 = np.pi*5**2

# Estágio 2
Isp_2 = 420  # s
F_2 = 4450  # kN
m0_2 = 677000  # kg
mp_2 = 429000  # kg
me_2 = 32600  # kg
epsilon_2 = 0.071
lambda_2 = 0.466
burn_time_2 = 360  # s
area_cross_2 = np.pi*5**2

# Estágio 3
Isp_3 = 420  # s
F_3 = 890  # kN
m0_3 = 215000  # kg
mp_3 = 109000  # kg
me_3 = 25700  # kg
epsilon_3 = 0.191
lambda_3 = 0.603
burn_time_3 = np.array([165, 335])  # s
area_cross_3 = np.pi*3.3**2

# Vetores numpy
Isp = np.array([Isp_1, Isp_2, Isp_3])
F = np.array([F_1, F_2, F_3])
m0 = np.array([m0_1, m0_2, m0_3])
mp = np.array([mp_1, mp_2, mp_3])
me = np.array([me_1, me_2, me_3])
epsilon = np.array([epsilon_1, epsilon_2, epsilon_3])
lambda_ = np.array([lambda_1, lambda_2, lambda_3])
area_cross = np.array([area_cross_1, area_cross_2, area_cross_3])
burn_time = np.array([burn_time_1, burn_time_2, burn_time_3[0], burn_time_3[1]])

# Exibindo os vetores
print("Impulso específico (s):", Isp)
print("Empuxo total (kN):", F)
print("Massa inicial total (kg):", m0)
print("Massa de propelente (kg):", mp)
print("Massa da estrutura e motores (kg):", me)
print("Razão estrutural:", epsilon)
print("Razão de carga útil:", lambda_)
print("Área da seção transversal (m^2):", area_cross)
print("Tempo de queima (s):", burn_time)


json_modelo_atm_terrestre = r'C:\Users\gt_po\Documents\tcc\mecvooespacial\src\domain\modelos\planeta\dados_JSON_planetas\dados_atmosfericos_terra.json'
terra = (ConstrutorPlaneta()
         .com_delta_temperatura_atm(10)
         .com_raio_equatorial(6378.1370e3)
         .com_velocidade_inercial_de_rotacao(7.2921150e-5)
         .com_gravidade(9.80665)
         .com_mut(3.986004418e14)
         .com_J2(0.00108263)
         .com_J3(-0.00000254)
         .com_J4(-0.00000161)
         .com_tempo_longitude_celeste_nula(0)
         .com_modelo_atmosferico(ModeloAtmosferico(json_modelo_atm_terrestre))
         .construir())

meuModeloPropulsivo = ConstrutorModeloPropulsivo().com_impulso_especifico(
    [290, 420, 420]).com_massa_propelente_estagios(
    [1997000,  429000,  109000]).com_duracao_queima_estagios(
    [168, 360, 165+335]).com_tempo_primeira_queima_terceiro_estagio(165).com_tempo_espera_separacao(
    [2, 2, 2]).com_tempo_espera_ignicao([2, 10]).com_massa_de_carga_util(81000).com_h0(0.0).com_planeta(terra).construir()

meuModeloEstrutural = ConstrutorModeloEstrutural().com_massa_estrutural_por_estagio(
    [137000 ,  40100,  15200 ]).com_massa_de_carga_util(
    81000).com_area_secao_transversal_1_estagio(76).com_area_secao_transversal_2_estagio(
    76).com_area_secao_transversal_3_estagio(35).com_area_secao_transversal_carga_util(
    30).com_comprimento_carga_util(1).com_comprimento_total_do_foguete(
    42 + 24.8 + 18.8).com_comprimento_sem_1_estagio(24.8 + 18.8).com_comprimento_sem_2_estagio(
    18.8).construir()

fogueteConceitual = ConstrutorDeFoguete().com_modelo_propulsivo(meuModeloPropulsivo).com_modelo_estrutural(
    meuModeloEstrutural).com_modelo_aerodinamico(
    ModeloAerodinamico()).construir()


alcantara = ConstrutorBaseDeLancamento().com_altitude_base(
    0).com_latitude_inicial().com_longitude_inicial().com_comprimento_trilho(
    fogueteConceitual.modelo_estrutural.comprimento_total_do_foguete).construir()

fogueteConceitual.mostra_dados()

# Condições iniciais
tempo_simulacao = 5000
velocidade_inicial = 0.5
angulo_elevacao_inicial = np.deg2rad(85)

# Criar uma órbita alvo
orbita_alvo = Orbita.circular(42.164140e6, np.deg2rad(5))
condicoes_iniciais = [tempo_simulacao, velocidade_inicial, angulo_elevacao_inicial, orbita_alvo]
simulacao = Simulacao(terra, alcantara, fogueteConceitual, condicoes_iniciais)

with open('../../exemplos/exemplos_simulacao_voo/simulacao.pkl', 'wb') as f:
    pickle.dump(simulacao, f)
resposta = simulacao.simular()






