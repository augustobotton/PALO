# -*- coding: utf-8 -*-


## Aplicação prática de cálculo de impulso de velocidade total
# Dados das 4 configurações de foguetes apresentados na referência
# Pedro L. K. da Cás, Carlos A. G. Veras, Olexiy Shynkarenko, and Rodrigo Leonardi. A
# brazilian space launch system for the small satellite market. Aerospace, 6(11), 2019
# VLM-1 e 3 veículos feitos com combinações de seus motores
#  1° estágio  |  2º estágio  |  3° estágio
# VLM-1 1xS50  |    1xS50     |    1xS44
#  C-1: 3xS50  |    1xS50     |    1xS44
#  C-2: 3xS50  |    1xS50     |    1xRD843
#  C-3: 3xS50  |    1xS50     |    4xRD843
#
import numpy as np
#
## Constantes
g = 9.81  # m/s**2
## Dados
# Massa de carga útil - Supor um CubeSat 8U (8*1,33+2,36 kg)
mL = 13  # kg
# Massa estutural
ms_vlm = np.array([1367, 1367, 166.5])  # kg - VLM-1
ms_c1 = np.array([4650, 1367, 166.5])  # kg - C-1
ms_c2 = np.array([4650, 1367, 161.9])  # kg - C-2
ms_c3 = np.array([4650, 1367, 228.7])  # kg - C-3
# Massa de propelente
mp_vlm = np.array([11058, 11058, 813])  # kg - VLM-1
mp_c1 = np.array([33157, 11058, 813])  # kg - C-1
mp_c2 = np.array([33157, 11058, 609])  # kg - C-2
mp_c3 = np.array([33157, 11058, 811])  # kg - C-3
# Impulso específico ideal
Isp_vlm = np.array([271, 271, 270])  # s - VLM-1
Isp_c1 = np.array([251, 271, 270])  # s - C-1
Isp_c2 = np.array([251, 271, 315])  # s - C-2
Isp_c3 = np.array([251, 271, 315])  # s - C-3
## Cálculos
# Razões estruturais
sigma_vlm = ms_vlm/(ms_vlm+mp_vlm)
sigma_c1 = ms_c1/(ms_c1+mp_c1)
sigma_c2 = ms_c2/(ms_c2+mp_c2)
sigma_c3 = ms_c3/(ms_c3+mp_c3)
# Massa total no inicio da queima de cada estagio
m0_vlm = np.empty(3); m0_c1 = np.empty(3); m0_c2 = np.empty(3); m0_c3 = np.empty(3)
# VLM-1
m0_vlm[0] = np.sum(ms_vlm)+np.sum(mp_vlm)+mL
m0_vlm[1] = m0_vlm[0]-ms_vlm[0]-mp_vlm[0]
m0_vlm[2] = m0_vlm[1]-ms_vlm[1]-mp_vlm[1]
# C-1
m0_c1[0] = np.sum(ms_c1)+np.sum(mp_c1)+mL
m0_c1[1] = m0_c1[0]-ms_c1[0]-mp_c1[0]
m0_c1[2] = m0_c1[1]-ms_c1[1]-mp_c1[1]
# C-2
m0_c2[0] = np.sum(ms_c2)+np.sum(mp_c2)+mL
m0_c2[1] = m0_c2[0]-ms_c2[0]-mp_c2[0]
m0_c2[2] = m0_c2[1]-ms_c2[1]-mp_c2[1]
# C-3
m0_c3[0] = np.sum(ms_c3)+np.sum(mp_c3)+mL
m0_c3[1] = m0_c3[0]-ms_c3[0]-mp_c3[0]
m0_c3[2] = m0_c3[1]-ms_c3[1]-mp_c3[1]
#
# Razões de carga útil
lamb_vlm = np.empty(3); lamb_c1 = np.empty(3); lamb_c2 = np.empty(3); lamb_c3 = np.empty(3)
# VLM-1
lamb_vlm[0] = m0_vlm[1]/m0_vlm[0]
lamb_vlm[1] = m0_vlm[2]/m0_vlm[1]
lamb_vlm[2] = mL/m0_vlm[2]
# C-1
lamb_c1[0] = m0_c1[1]/m0_c1[0]
lamb_c1[1] = m0_c1[2]/m0_c1[1]
lamb_c1[2] = mL/m0_c1[2]
# C-2
lamb_c2[0] = m0_c2[1]/m0_c2[0]
lamb_c2[1] = m0_c2[2]/m0_c2[1]
lamb_c2[2] = mL/m0_c2[2]
# C-3
lamb_c3[0] = m0_c3[1]/m0_c3[0]
lamb_c3[1] = m0_c3[2]/m0_c3[1]
lamb_c3[2] = mL/m0_c3[2]
# Razão de carga útil total
# VLM-1
lambL_vlm = np.prod(lamb_vlm)
# C-1
lambL_c1 = np.prod(lamb_c1)
# C-2
lambL_c2 = np.prod(lamb_c2)
# C-3
lambL_c3 = np.prod(lamb_c3)
# Impulso de velocidade
# VLM-1
Dv_vlm = -np.sum(g*Isp_vlm*np.log(sigma_vlm+(1-sigma_vlm)*lamb_vlm))
# C-1
Dv_c1 = -np.sum(g*Isp_c1*np.log(sigma_c1+(1-sigma_c1)*lamb_c1))
# C-2
Dv_c2 = -np.sum(g*Isp_c2*np.log(sigma_c2+(1-sigma_c2)*lamb_c2))
# C-3
Dv_c3 = -np.sum(g*Isp_c3*np.log(sigma_c3+(1-sigma_c3)*lamb_c3))
## Resultados
print('*********************************************')
print('VLM-1')
print('*********************************************')
print('Massas iniciais antes da queima de cada estagio - kg'); print(m0_vlm)
print('Razoes estruturais'); print(sigma_vlm)
print('Razoes de carga util'); print(lamb_vlm)
print('Razao de carga útil total'); print(lambL_vlm)
print('Impulso de velocidade total - km/s'); print(Dv_vlm/1e3)
print('*********************************************')
print('C-1')
print('*********************************************')
print('Massas iniciais antes da queima de cada estagio - kg'); print(m0_c1)
print('Razoes estruturais'); print(sigma_c1)
print('Razoes de carga util'); print(lamb_c1)
print('Razao de carga útil total'); print(lambL_c1)
print('Impulso de velocidade total - km/s'); print(Dv_c1/1e3)
print('*********************************************')
print('C-2')
print('*********************************************')
print('Massas iniciais antes da queima de cada estagio - kg'); print(m0_c2)
print('Razoes estruturais'); print(sigma_c2)
print('Razoes de carga util'); print(lamb_c2)
print('Razao de carga útil total'); print(lambL_c2)
print('Impulso de velocidade total - km/s'); print(Dv_c2/1e3)
print('*********************************************')
print('C-3')
print('*********************************************')
print('Massas iniciais antes da queima de cada estagio - kg'); print(m0_c3)
print('Razoes estruturais'); print(sigma_c3)
print('Razoes de carga util'); print(lamb_c3)
print('Razao de carga útil total'); print(lambL_c3)
print('Impulso de velocidade total - km/s'); print(Dv_c3/1e3)
