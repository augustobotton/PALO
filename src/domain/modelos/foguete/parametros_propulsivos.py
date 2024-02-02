# Parâmetros propulsivos
import numpy as np

from src.domain.request import parametros

impulso_especico_por_estagio = np.array([251, 271, 315])  # s - Impulso específico dos estágios
parametros.impulso_especifico_por_estagio = impulso_especico_por_estagio
mp[0] = 55290
parametros.mp[0] = mp[0]
mp[1] = 11058  # kg - Massa de propelente dos estágios
parametros.mp[1] = mp[1]
parametros.mp = mp
mp3 = 224.53  # kg - Massa de propelente do terceiro estágio (RD843)
parametros.mp3 = mp3
Tq1 = 62.0233  # s - DADO DOS MOTORES DO 1° ESTÁGIO
parametros.Tq1 = Tq1
Tq2 = 64.6105  # s - DADO DOS MOTORES DO 2° ESTÁGIO
parametros.Tq2 = Tq2
Tq3 = 277.5325  # s - TEMPO DE QUEIMA DO 3° ESTÁGIO SE ELE IGNITASSE SÓ UMA VEZ
parametros.Tq3 = Tq3

# Parâmetros de massa estrutural e de carga útil
massa_estrutural_por_estagio = np.array([7385, 1367, 59.69])  # kg - Massa estrutural dos estágios
parametros.massa_estrutural_por_estagio = massa_estrutural_por_estagio
massa_de_carga_util = 13  # kg - Massa da carga útil
parametros.massa_de_carga_util = massa_de_carga_util

# Parâmetros aerodinâmicos e ambientais
fator_correcao = 1.28  # Fator de correção do arrasto
parametros.fator_correcao_arrasto = fator_correcao
S1 = 4.6 * 5 / 3  # m^2 - Area aproximada da secao transversal do primeiro estagio
S2 = 1.5  # m^2 - Area aproximada da secao transversal do segundo estagio
S3 = 1.5  # m^2 - Area aproximada da secao transversal do terceiro estagio
SL = 1.5  # m^2 - Area aproximada da secao transversal da carga útil

lt = 7.33 + 7.1 + 6.28  # m - Comprimento total
l2 = 7.1 + 6.28  # Comprimento sem o primeiro estagio
l3 = 6.28  # Comprimento sem o segundo estagio

l4 = 1  # Comprimento da carga util
f2 = (l2 / lt) * 0.5 + 0.5  # Fator de correcao do segundo estagio
f3 = (l3 / lt) * 0.5 + 0.5  # Fator de correcao do terceiro estagio
f4 = (l4 / lt) * 0.5 + 0.5  # Fator de correcao da carga util

Sr = np.array([S1, S2 * f2, S3 * f3, SL * f4])
Sr = Sr.reshape(-1, 1)
parametros.area_de_referencia = Sr

lc = 1.5  # Comprimento característico - diâmetro dos estágios 2 e superiores
parametros.lc = lc
