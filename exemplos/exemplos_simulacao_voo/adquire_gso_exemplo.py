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
    [251, 271, 315]).com_massa_propelente_estagios(
    [77366, 11058, 225.33]).com_duracao_queima_estagios(
    [62, 64.62, 279]).com_tempo_primeira_queima_terceiro_estagio(213).com_tempo_espera_separacao(
    [2, 2, 2]).com_tempo_espera_ignicao([5, 580]).com_massa_de_carga_util(13).com_h0(0.0).com_planeta(terra).construir()

meuModeloEstrutural = ConstrutorModeloEstrutural().com_massa_estrutural_por_estagio(
    [7750, 1367, 64.7544]).com_massa_de_carga_util(
    13).com_area_secao_transversal_1_estagio(4.6 * 5 / 3).com_area_secao_transversal_2_estagio(
    1.5).com_area_secao_transversal_3_estagio(1.5).com_area_secao_transversal_carga_util(
    1.5).com_comprimento_carga_util(1).com_comprimento_total_do_foguete(
    7.33 + 7.1 + 6.28).com_comprimento_sem_1_estagio(7.1 + 6.28).com_comprimento_sem_2_estagio(
    6.28).construir()

fogueteConceitual = ConstrutorDeFoguete().com_modelo_propulsivo(meuModeloPropulsivo).com_modelo_estrutural(
    meuModeloEstrutural).com_modelo_aerodinamico(
    ModeloAerodinamico()).construir()


alcantara = ConstrutorBaseDeLancamento().com_altitude_base(
    0).com_latitude_inicial().com_longitude_inicial().com_comprimento_trilho(
    fogueteConceitual.modelo_estrutural.comprimento_total_do_foguete).construir()

fogueteConceitual.mostra_dados()

# Condições iniciais
tempo_simulacao = 120000
velocidade_inicial = 1
angulo_elevacao_inicial = np.deg2rad(80)

# Criar uma órbita alvo
orbita_alvo = Orbita.circular(42.164140e6, np.deg2rad(8))
condicoes_iniciais = [tempo_simulacao, velocidade_inicial, angulo_elevacao_inicial, orbita_alvo]
simulacao = Simulacao(terra, alcantara, fogueteConceitual, condicoes_iniciais)

with open('../../exemplos/exemplos_simulacao_voo/simulacao.pkl', 'wb') as f:
    pickle.dump(simulacao, f)
resposta = simulacao.simular()






