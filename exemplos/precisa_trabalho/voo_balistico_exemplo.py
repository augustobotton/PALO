import pickle

import numpy as np

from src.domain.modelos.Simulacao import Simulacao
from src.domain.modelos.foguete.aerodinamica.ModeloAerodinamico import ModeloAerodinamico
from src.domain.modelos.foguete.estrutura.ConstrutorModeloEstrutural import ConstrutorModeloEstrutural
from src.domain.modelos.foguete.propulsao.ConstrutorModeloPropulsivo import ConstrutorModeloPropulsivo
from src.domain.modelos.foguete.veiculo_lancador.ConstrutorFoguete import ConstrutorDeFoguete
from src.domain.modelos.orbitas import Orbita
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


meuModeloPropulsivo = ConstrutorModeloPropulsivo().com_impulso_especifico([260.6, 261.1]).com_massa_propelente_estagios(
    [677, 898]).com_duracao_queima_estagios([62, 64.62]).com_tempo_espera_separacao(
    [2, 2]).com_tempo_espera_ignicao([5]).com_massa_estrutural_por_estagio(
    [284, 320]).com_massa_de_carga_util(400).com_h0(0.0).com_planeta(terra).construir()

#TODO estou criando dois atributos para dois modelos diferentes, mas o que eu quero é criar um atributo que seja comum para os dois modelos.

meuModeloEstrutural = ConstrutorModeloEstrutural().com_massa_estrutural_por_estagio(
    [284, 320]).com_massa_de_carga_util(400).com_area_secao_transversal_1_estagio(np.pi*(0.557/2)**2).com_area_secao_transversal_2_estagio(
    np.pi*(0.557/2)**2).com_area_secao_transversal_carga_util(
    np.pi*(0.46/2)**2).com_comprimento_total_do_foguete(12.6).com_comprimento_sem_1_estagio(12.6-3.214).com_comprimento_carga_util(12.6-3.294).construir()

fogueteConceitual = ConstrutorDeFoguete().com_modelo_propulsivo(meuModeloPropulsivo).com_modelo_estrutural(
    meuModeloEstrutural).com_modelo_aerodinamico(
    ModeloAerodinamico()).construir()


alcantara = ConstrutorBaseDeLancamento().com_altitude_base(
    0).com_latitude_inicial().com_longitude_inicial().com_comprimento_trilho(
    fogueteConceitual.modelo_estrutural.comprimento_total_do_foguete).construir()

fogueteConceitual.mostra_dados()

# Condições iniciais
tempo_simulacao = 75000
velocidade_inicial = 1
angulo_elevacao_inicial = np.deg2rad(80)
phi_inicial = np.deg2rad(80)

# Criar uma órbita alvo
orbita_alvo = Orbita.circular(42.164140e6, np.deg2rad(5))
condicoes_iniciais = [tempo_simulacao, velocidade_inicial, angulo_elevacao_inicial, orbita_alvo, phi_inicial]
simulacao = Simulacao(terra, alcantara, fogueteConceitual, condicoes_iniciais)

with open('../../src/domain/construtorderesultados/simulacao.pkl', 'wb') as f:
    pickle.dump(simulacao, f)
resposta = simulacao.simular()
