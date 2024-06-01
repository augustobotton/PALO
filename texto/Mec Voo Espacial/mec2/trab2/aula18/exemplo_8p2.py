# Exemplo 8.2 do livro TEWARI, A. Atmospheric and Space Flight Dynamics: Modelling
# and simulation with MATLAB and Simulink. Boston: Birkhauser, 2007.
# Variaveis globais
# Para passagem de parametros

from scipy.optimize import root
import numpy as np

# Variáveis globais
Isp = np.array([290, 290, 455])  # s
sig = np.array([0.07, 0.07, 0.07])  # Razões estruturais
alf = np.array([1, 1.2, 0.65])  # Razão de carga útil de cada estágio normalizada pela do primeiro
mL = 1000  # kg - Massa de carga útil
Dv = 13000  # m/s - Impulso de velocidade total do foguete
g = 9.81  # m/s^2

# Função objetivo para encontrar a razão de carga útil do primeiro estágio
def obj_eq_fog(lam1):
    sdv = np.sum(-g * Isp * np.log(sig + (1 - sig) * lam1 * alf))
    y = Dv - sdv
    return y

# Determinação da razão de carga útil do primeiro estágio
lam1_chute_inicial = 0.5
sol = root(obj_eq_fog, lam1_chute_inicial)
lam1 = sol.x[0]

# Determinação da massa no início da queima de cada estágio
lam = lam1 * alf  # Razão de carga útil de cada estágio
lamT = lam.prod()  # Razão de carga útil total
m01 = mL / lamT
m02 = m01 * lam[0]
m03 = m02 * lam[1]

# Determinação da massa de propelente de cada estágio
msp = np.array([m01 - m02, m02 - m03, m03 - mL])  # Massa estrutural e total de propelente em cada estágio
ms = msp * sig  # Massa estrutural em cada estágio
mp = msp - ms  # Massa de propelente em cada estágio
mpT = mp.sum()  # Massa de propelente total
pmp = 100 * mpT / m01  # Percentual de massa de propelente com respeito à massa total do foguete
msT = ms.sum()  # Massa estrutural total
pms = 100 * msT / m01  # Percentual de massa estrutural com respeito à massa total do foguete

# Saída de resultados
print("Razão de carga útil de cada estágio: ")
print(f"lambda_1: {lam[0]}, lambda_2: {lam[1]}, lambda_3: {lam[2]}")
print(f"Razão de carga útil total: {lamT}")
print("Massa no início da queima de cada estágio: (kg)")
print(f"m01: {m01}, m02: {m02}, m03: {m03}")
print("Massa estrutural e total de propelente em cada estágio (kg): ")
print(f"Primeiro estágio: {msp[0]}, Segundo estágio: {msp[1]}, Terceiro estágio: {msp[2]}")
print("Massa de propelente em cada estágio (kg): ")
print(f"Primeiro estágio: {mp[0]}, Segundo estágio: {mp[1]}, Terceiro estágio: {mp[2]}")
print(f"Massa total de propelente (kg): {mpT}")
print(f"Percentual de massa de propelente com respeito à massa total do foguete: {pmp}")
print("Massa estrutural em cada estágio (kg): ")
print(f"Primeiro estágio: {ms[0]}, Segundo estágio: {ms[1]}, Terceiro estágio: {ms[2]}")
print(f"Massa estrutural total (kg): {msT}")
print(f"Percentual de massa estrutural com respeito à massa total do foguete: {pms}")
