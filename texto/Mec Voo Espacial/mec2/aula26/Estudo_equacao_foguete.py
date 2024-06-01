import numpy as np

# Constantes
g = 9.81  # m/s^2

# Dados
# Massa de carga útil - Supor um CubeSat 8U (8*1,33+2,36 kg)
mL = 13  # kg

# Massa estrutural
ms_vlm = np.array([1367, 1367, 166.5])  # kg - VLM-1
ms_c1 = np.array([4650*6/3, 1367, 166.5])  # kg - C-1
ms_c2 = np.array([4650*6/3, 1367, 161.9])  # kg - C-2
ms_c3 = np.array([4650*6/3, 1367, 228.7])  # kg - C-3

# Massa de propelente
mp_vlm = np.array([11058, 11058, 813])  # kg - VLM-1
mp_c1 = np.array([33157*6/3, 11058, 813])  # kg - C-1
mp_c2 = np.array([33157*6/3, 11058, 609])  # kg - C-2
mp_c3 = np.array([33157*6/3, 11058, 811])  # kg - C-3

# Impulso específico ideal do primeiro e segundo estágio
Isp_vlm = np.array([271, 271, 270])  # s - VLM-1
Isp_c1 = np.array([251, 271, 270])  # s - C-1
Isp_c2 = np.array([251, 271, 315])  # s - C-2
Isp_c3 = np.array([251, 271, 315])  # s - C-3

# Cálculos
# Razões estruturais
sigma_vlm = ms_vlm / (ms_vlm + mp_vlm)
sigma_c1 = ms_c1 / (ms_c1 + mp_c1)
sigma_c2 = ms_c2 / (ms_c2 + mp_c2)
sigma_c3 = ms_c3 / (ms_c3 + mp_c3)

# Funções de massa inicial e razão de carga útil
def calc_m0(ms, mp):
    m0 = np.zeros(3)
    m0[0] = sum(ms) + sum(mp) + mL
    m0[1] = m0[0] - ms[0] - mp[0]
    m0[2] = m0[1] - ms[1] - mp[1]
    return m0

def calc_lamb(m0):
    lamb = np.zeros(3)
    lamb[0] = m0[1] / m0[0]
    lamb[1] = m0[2] / m0[1]
    lamb[2] = mL / m0[2]
    return lamb

m0_vlm = calc_m0(ms_vlm, mp_vlm)
m0_c1 = calc_m0(ms_c1, mp_c1)
m0_c2 = calc_m0(ms_c2, mp_c2)
m0_c3 = calc_m0(ms_c3, mp_c3)

lamb_vlm = calc_lamb(m0_vlm)
lamb_c1 = calc_lamb(m0_c1)
lamb_c2 = calc_lamb(m0_c2)
lamb_c3 = calc_lamb(m0_c3)

# Razão de carga útil total
lambL_vlm = np.prod(lamb_vlm)
lambL_c1 = np.prod(lamb_c1)
lambL_c2 = np.prod(lamb_c2)
lambL_c3 = np.prod(lamb_c3)

# Impulso de velocidade
def calc_Dv(Isp, sigma, lamb):
    return -np.sum(g * Isp * np.log(sigma + (1 - sigma) * lamb))

Dv_vlm = calc_Dv(Isp_vlm, sigma_vlm, lamb_vlm)
Dv_c1 = calc_Dv(Isp_c1, sigma_c1, lamb_c1)
Dv_c2 = calc_Dv(Isp_c2, sigma_c2, lamb_c2)
Dv_c3 = calc_Dv(Isp_c3, sigma_c3, lamb_c3)

# Resultados
print('*********************************************')
print('VLM-1')
print('*********************************************')
print('Massas iniciais antes da queima de cada estagio - kg');print(m0_vlm)
print('Massa total do primeiro estagio (kg): ', m0_vlm[0] - m0_vlm[1])
print('Massa total do segundo estagio (kg): ', m0_vlm[1] - m0_vlm[2])
print('Massa total do terceiro estagio (kg): ', m0_vlm[2] - mL)
print('Razoes estruturais', sigma_vlm)
print('Razoes de carga util', lamb_vlm)
print('Razao de carga útil total', lambL_vlm)
print('Impulso de velocidade total - km/s', Dv_vlm / 1e3)
# Resultados
print('*********************************************')
print('C1')
print('*********************************************')
print('Massas iniciais antes da queima de cada estagio - kg');print(m0_c1)
print('Massa total do primeiro estagio (kg): ', m0_c1[0] - m0_c1[1])
print('Massa total do segundo estagio (kg): ', m0_c1[1] - m0_c1[2])
print('Massa total do terceiro estagio (kg): ', m0_c1[2] - mL)
print('Razoes estruturais', sigma_c1)
print('Razoes de carga util', lamb_c1)
print('Razao de carga útil total', lambL_c1)
print('Impulso de velocidade total - km/s', Dv_c1 / 1e3)
# Resultados
print('*********************************************')
print('C2')
print('*********************************************')
print('Massas iniciais antes da queima de cada estagio - kg');print(m0_c2)
print('Massa total do primeiro estagio (kg): ', m0_c2[0] - m0_c2[1])
print('Massa total do segundo estagio (kg): ', m0_c2[1] - m0_c2[2])
print('Massa total do terceiro estagio (kg): ', m0_c2[2] - mL)
print('Razoes estruturais', sigma_c2)
print('Razoes de carga util', lamb_c2)
print('Razao de carga útil total', lambL_c2)
print('Impulso de velocidade total - km/s', Dv_c2 / 1e3)
# Resultados
print('*********************************************')
print('VLM-1')
print('*********************************************')
print('Massas iniciais antes da queima de cada estagio - kg');print(m0_c3)
print('Massa total do primeiro estagio (kg): ', m0_c3[0] - m0_c3[1])
print('Massa total do segundo estagio (kg): ', m0_c3[1] - m0_c3[2])
print('Massa total do terceiro estagio (kg): ', m0_c3[2] - mL)
print('Razoes estruturais', sigma_c3)
print('Razoes de carga util', lamb_c3)
print('Razao de carga útil total', lambL_c3)
print('Impulso de velocidade total - km/s', Dv_c3 / 1e3)

