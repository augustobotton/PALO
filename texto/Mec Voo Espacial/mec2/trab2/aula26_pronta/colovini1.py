import numpy as np
import matplotlib.pyplot as plt

## Aplicação prática de cálculo de impulso de velocidade total
# Dados apresentados na referência
# Pedro L. K. da Cás, Carlos A. G. Veras, Olexiy Shynkarenko, and Rodrigo Leonardi. A
# brazilian space launch system for the small satellite market. Aerospace, 6(11), 2019
# Veiculo  C-2  da referencia alterado no terceiro estagio. Chutar
# massas de propelente. Configuracao com 3, 4 e 5 motores de primeiro estagio.

#       1° estágio| 2º estágio | 3° estágio
# C-2-c:  3xS50     | 1xS50      | 1xRD843
# C-2-d:  4xS50     | 1xS50      | 1xRD843
# C-2-e:  5xS50     | 1xS50      | 1xRD843


# Constantes
g = 9.81  # m/s^2

# Dados
# Massa de carga útil - Supor um CubeSat 8U (8*1,33+2,36 kg)
mL = 13  # kg

# Variacoes da massa de propelente
fp = np.arange(0.1, 2.01, 0.01)  # De metade até o dobro
N = len(fp)
Dv_c2c = np.zeros(N)
Dv_c2d = np.zeros(N)
Dv_c2e = np.zeros(N)  # Impulsos resultantes
ms3 = np.zeros(N)  # Massa estrutural do 3° estagio

for i in range(N):
    # Massa de propelente
    mp_c2c = np.array([33157, 11058, fp[i] * 609])  # kg - C-2-c
    mp_c2d = np.array([33157 * 4 / 3, 11058, fp[i] * 609])  # kg - C-2-d
    mp_c2e = np.array([33157 * 5 / 3, 11058, fp[i] * 609])  # kg - C-2-e

    # Massa estrutural
    ms3[i] = (0.21 / (1 - 0.21)) * mp_c2c[2]
    ms_c2c = np.array([4650, 1367, ms3[i]])  # kg - C-2-c
    ms_c2d = np.array([4650 * 4 / 3, 1367, ms3[i]])  # kg - C-2-d
    ms_c2e = np.array([4650 * 5 / 3, 1367, ms3[i]])  # kg - C-2-e

    # Impulso específico ideal do primeiro e segundo estágio
    Isp_c2 = np.array([251, 271, 315])  # s - C-2

    # Cálculos
    # Razões estruturais
    sigma_c2c = ms_c2c / (ms_c2c + mp_c2c)
    sigma_c2d = ms_c2d / (ms_c2d + mp_c2d)
    sigma_c2e = ms_c2e / (ms_c2e + mp_c2e)

    # Massa total no início da queima de cada estágio
    # C-2-c
    m0_c2c = np.zeros(3)
    m0_c2c[0] = np.sum(ms_c2c) + np.sum(mp_c2c) + mL
    m0_c2c[1] = m0_c2c[0] - ms_c2c[0] - mp_c2c[0]
    m0_c2c[2] = m0_c2c[1] - ms_c2c[1] - mp_c2c[1]

    # C-2-d
    m0_c2d = np.zeros(3)
    m0_c2d[0] = np.sum(ms_c2d) + np.sum(mp_c2d) + mL
    m0_c2d[1] = m0_c2d[0] - ms_c2d[0] - mp_c2d[0]
    m0_c2d[2] = m0_c2d[1] - ms_c2d[1] - mp_c2d[1]

    # C-2-e
    m0_c2e = np.zeros(3)
    m0_c2e[0] = np.sum(ms_c2e) + np.sum(mp_c2e) + mL
    m0_c2e[1] = m0_c2e[0] - ms_c2e[0] - mp_c2e[0]
    m0_c2e[2] = m0_c2e[1] - ms_c2e[1] - mp_c2e[1]

    # Razões de carga útil
    # C-2-c
    lamb_c2c = np.zeros(3)
    lamb_c2c[0] = m0_c2c[1] / m0_c2c[0]
    lamb_c2c[1] = m0_c2c[2] / m0_c2c[1]
    lamb_c2c[2] = mL / m0_c2c[2]

    # C-2-d
    lamb_c2d = np.zeros(3)
    lamb_c2d[0] = m0_c2d[1] / m0_c2d[0]
    lamb_c2d[1] = m0_c2d[2] / m0_c2d[1]
    lamb_c2d[2] = mL / m0_c2d[2]

    # C-2-e
    lamb_c2e = np.zeros(3)
    lamb_c2e[0] = m0_c2e[1] / m0_c2e[0]
    lamb_c2e[1] = m0_c2e[2] / m0_c2e[1]
    lamb_c2e[2] = mL / m0_c2e[2]

    # Razão de carga útil total
    # C-2-c
    lambL_c2c = np.prod(lamb_c2c)
    # C-2-d
    lambL_c2d = np.prod(lamb_c2d)
    # C-2-e
    lambL_c2e = np.prod(lamb_c2e)

    # Impulso de velocidade
    # C-2-c
    Dv_c2c[i] = -np.sum(g * Isp_c2 * np.log(sigma_c2c + (1 - sigma_c2c) * lamb_c2c))
    # C-2-d
    Dv_c2d[i] = -np.sum(g * Isp_c2 * np.log(sigma_c2d + (1 - sigma_c2d) * lamb_c2d))
    # C-2-e
    Dv_c2e[i] = -np.sum(g * Isp_c2 * np.log(sigma_c2e + (1 - sigma_c2e) * lamb_c2e))

# Plotting


fig, axes = plt.subplots(1, 2, figsize=(8, 8))

# Plot the first three subplots in the first axis
axes[0].plot(fp * 100, Dv_c2c / 1e3)
axes[0].plot(fp * 100, Dv_c2d / 1e3)
axes[0].plot(fp * 100, Dv_c2e / 1e3)
axes[0].set_ylabel('Δv - km/s')
axes[0].set_xlabel('Fração de m_p original - %')
axes[0].grid()
axes[0].legend(['3 motores S50', '4 motores S50', '5 motores S50'])

# Plot the fourth subplot in the second axis
axes[1].plot(fp * 100, ms3)
axes[1].set_xlabel('Fração de m_p original - %')
axes[1].set_ylabel('m_{s_3} - kg')
axes[1].grid()

plt.subplots_adjust(hspace=0.25)  # Adjust spacing between subplots

plt.show()