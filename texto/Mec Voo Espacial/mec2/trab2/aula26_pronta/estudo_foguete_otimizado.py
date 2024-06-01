import numpy as np
import matplotlib.pyplot as plt

# Constantes
g = 9.81  # m/s^2

# Massa de carga útil - Supor um CubeSat 8U (8*1,33+2,36 kg)
mL = 13  # kg

# Variacoes da massa de propelente
fp = np.arange(0.1, 2.01, 0.01)  # De metade ate o dobro
N = len(fp)
Dv_c2c = np.zeros(N)
Dv_c2d = np.zeros(N)
Dv_c2e = np.zeros(N)
Dv_c2f = np.zeros(N)
ms3 = np.zeros(N)  # Massa estrutural do 3° estagio
ms4 = np.zeros(N)  # Massa estrutural do 7 motores

for i in range(N):
    # Massa de propelente
    mp_c2c = [33157, 11058, fp[i]*609]  # kg - C-2-c
    mp_c2d = [33157*4/3, 11058, fp[i]*609]  # kg - C-2-d
    mp_c2e = [33157*5/3, 11058, fp[i]*609]  # kg - C-2-e
    mp_c2f = [33157*7/3, 11058, fp[i]*609]  # kg - C-2-f

    # Massa estutural
    ms3[i] = (0.21/(1-0.21))*mp_c2c[2]
    ms4[i] = (0.21/(1-0.21))*mp_c2f[2]
    ms_c2c = [4650, 1367, ms3[i]]  # kg - C-2-c
    ms_c2d = [4650*4/3, 1367, ms3[i]]  # kg - C-2-d
    ms_c2e = [4650*5/3, 1367, ms3[i]]  # kg - C-2-e
    ms_c2f = [4650*7/3, 1367, ms4[i]]  # kg - C-2-f

    # Impulso específico ideal do primeiro e segundo estágio
    Isp_c2 = [251, 271, 315]  # s - C-2

    # Razões estruturais
    sigma_c2c = np.divide(ms_c2c, np.add(ms_c2c, mp_c2c))
    sigma_c2d = np.divide(ms_c2d, np.add(ms_c2d, mp_c2d))
    sigma_c2e = np.divide(ms_c2e, np.add(ms_c2e, mp_c2e))
    sigma_c2f = np.divide(ms_c2f, np.add(ms_c2f, mp_c2f))

    # Massa total no inicio da queima de cada estagio
    m0_c2c = [sum(ms_c2c)+sum(mp_c2c)+mL,
              sum(ms_c2c[1:])+sum(mp_c2c[1:])+mL,
              ms_c2c[2]+mp_c2c[2]+mL]

    m0_c2d = [sum(ms_c2d)+sum(mp_c2d)+mL,
              sum(ms_c2d[1:])+sum(mp_c2d[1:])+mL,
              ms_c2d[2]+mp_c2d[2]+mL]

    m0_c2e = [sum(ms_c2e)+sum(mp_c2e)+mL,
              sum(ms_c2e[1:])+sum(mp_c2e[1:])+mL,
              ms_c2e[2]+mp_c2e[2]+mL]
              
    m0_c2f = [sum(ms_c2f)+sum(mp_c2f)+mL,
              sum(ms_c2f[1:])+sum(mp_c2f[1:])+mL,
              ms_c2f[2]+mp_c2f[2]+mL]

    # Razões de carga útil
    lamb_c2c = [m0_c2c[1]/m0_c2c[0], m0_c2c[2]/m0_c2c[1], mL/m0_c2c[2]]
    lamb_c2d = [m0_c2d[1]/m0_c2d[0], m0_c2d[2]/m0_c2d[1], mL/m0_c2d[2]]
    lamb_c2e = [m0_c2e[1]/m0_c2e[0], m0_c2e[2]/m0_c2e[1], mL/m0_c2e[2]]
    lamb_c2f = [m0_c2f[1]/m0_c2f[0], m0_c2f[2]/m0_c2f[1], mL/m0_c2f[2]]

    # Razão de carga útil total
    lambL_c2c = np.prod(lamb_c2c)
    lambL_c2d = np.prod(lamb_c2d)
    lambL_c2e = np.prod(lamb_c2e)
    lambL_c2f = np.prod(lamb_c2f)
    
    # Impulso de velocidade
    Dv_c2c[i] = -sum(g*np.array(Isp_c2)*np.log(sigma_c2c+(1-sigma_c2c)*np.array(lamb_c2c)))
    Dv_c2d[i] = -sum(g*np.array(Isp_c2)*np.log(sigma_c2d+(1-sigma_c2d)*np.array(lamb_c2d)))
    Dv_c2e[i] = -sum(g*np.array(Isp_c2)*np.log(sigma_c2e+(1-sigma_c2e)*np.array(lamb_c2e)))
    Dv_c2f[i] = -sum(g*np.array(Isp_c2)*np.log(sigma_c2f+(1-sigma_c2f)*np.array(lamb_c2f)))

    # Cálculo da Tração
    tracao_c2c = (Isp_c2[0] * mp_c2c[0] * g) / 1000
    tracao_c2d = (Isp_c2[0] * mp_c2d[0] * g) / 1000
    tracao_c2e = (Isp_c2[0] * mp_c2e[0] * g) / 1000
    tracao_c2f = (Isp_c2[0] * mp_c2f[0] * g) / 1000

    

  # Print results
    print(f"Configuração C-2-c:\nMassa de Propelente: {mp_c2c}\nMassa Estrutural: {ms_c2c}\nImpulso Específico (ISP): {Isp_c2}\nTração: {tracao_c2c} kN\n")
    print(f"Configuração C-2-d:\nMassa de Propelente: {mp_c2d}\nMassa Estrutural: {ms_c2d}\nImpulso Específico (ISP): {Isp_c2}\nTração: {tracao_c2d} kN\n")
    print(f"Configuração C-2-e:\nMassa de Propelente: {mp_c2e}\nMassa Estrutural: {ms_c2e}\nImpulso Específico (ISP): {Isp_c2}\nTração: {tracao_c2e} kN\n")
    print(f"Configuração C-2-f:\nMassa de Propelente: {mp_c2f}\nMassa Estrutural: {ms_c2f}\nImpulso Específico (ISP): {Isp_c2}\nTração: {tracao_c2f} kN\n")

# Plots
fig, axs = plt.subplots(2, 1, figsize=(10, 12))

axs[0].plot(fp*100, Dv_c2c/1e3, label='3 motores S50')
axs[0].plot(fp*100, Dv_c2d/1e3, label='4 motores S50')
axs[0].plot(fp*100, Dv_c2e/1e3, label='5 motores S50')
axs[0].plot(fp*100, Dv_c2f/1e3, label='7 motores S50')
axs[0].set_ylabel('\Delta v - km/s')
axs[0].grid(True)
axs[0].legend()

axs[1].plot(fp*100, ms3, label='Massa estrutural do 3º estágio com 3 motores S50')
axs[1].plot(fp*100, ms4, label='Massa estrutural do 3º estágio com 7 motores S50')
axs[1].set_xlabel('Variação na massa de propelente do terceiro estágio - %')
axs[1].set_ylabel('Massa estrutural - kg')
axs[1].grid(True)
axs[1].legend()

plt.tight_layout()
plt.show()

