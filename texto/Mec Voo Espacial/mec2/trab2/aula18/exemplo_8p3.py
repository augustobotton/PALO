import numpy as np

# Dados do exemplo 8.3
mpb = 30000  
Ispb = 200   
sigb = 0.05  
mp10 = 25000 

# Resultados do exemplo 8.2
m01 = 77847.7296  
m02 = 19806.1060
m03 = 6046.9099
sigk = 0.07
sig = [0, sigk, sigk]

lam = [0, 0, 0.305305336, 0.165373724]

Isp1 = 290
Isp2 = 290
Isp3 = 455
mp3 = 4693.62617

# Cálculos
ms1 = sigk * (m01 - m02)
print('Massa estrutural do primeiro estagio do veiculo nucleo:')
print(ms1)

mp1 = m01 - m02 - ms1
print('Massa de propelente do primeiro estagio do veiculo nucleo:')
print(mp1)

msb = (sigb / (1 - sigb)) * mpb
print('Massa estrutural dos boosters:')
print(msb)

m00 = m01 + mpb + msb
print('Massa inicial do foguete com os boosters:')
print(m00)

sig[0] = (msb + ms1) / (msb + ms1 + mpb + mp10)
print('Razao estrutural do estagio zero:')
print(sig[0])

lam[0] = (m01 - mp10) / m00
print('Razao de carga util do estagio zero:')
print(lam[0])

sig[1] = ms1 / (ms1 + mp1 - mp10)
print('Razao estrutural do primeiro estagio modificado:')
print(sig[1])

lam[1] = m02 / (m01 - mp10)
print('Razao de carga util do primeiro estagio modificado:')
print(lam[1])

ve = [0, 0, 0, 0]
veb = 9.81*Ispb 
ve[1] = Isp1 * 9.81
ve[0] = (mpb * Ispb * 9.81 + mp10 * ve[1]) / (mpb + mp10)


print('Velocidade de exaustao media do estagio zero:')
print(ve[0])

ve[2]=Isp2*9.81
ve[3]=Isp3*9.81

Dv = -sum(ve[i] * np.log(sig[i] + (1 - sig[i]) * lam[i]) for i in range(3))
print('Impulso de velocidade total:')
print(Dv)

DDv = Dv - 13000
print('Acrescimo de impulso de velocidade (m/s):')
print(DDv)

des = 100 * DDv / 13000
print('Variacao percentual do desempenho:')
print(des)

lamT = np.prod(lam)
print('Nova razao de carga util total:')
print(lamT)

DlamT = lamT - 0.0128
print('Variacao da razao de carga util total:')
print(DlamT)

ef = 100 * DlamT / 0.0128
print('Variacao percentual da eficiencia:')
print(ef)

k = - (sum(ve[i] * np.log(sig[i] + (1 - sig[i]) * lam[i]) for i in range(3)) + 13000) / ve[3]
c1 = np.exp(k)
ms3 = m03 - mp3 - 1000
c2 = ms3 / m03

sig[2] = c2 / (1 - c1 + c2)
print('Nova razao estrutural do 3� estagio para Dv constante:')
print(sig[2])

lam[3] = c1 - c2
print('Nova razao de carga util do 3� estagio para Dv constante:')
print(lam[3])

mL = lam[3] * m03
print('Nova massa de carga util (kg):')
print(mL)

Dm = mL - 1000
print('Aumento de carga util (kg):')
print(Dm)

print('Aumento percentual de carga util:')
print(100 * (mL - 1000) / 1000)

mp3 = mp3 - Dm
print('Nova massa de propelente do terceiro estagio (kg):')
print(mp3)

lamTn = np.prod(lam)
print('Nova razao de carga util total:')
print(lamTn)

DlamT = lamTn - 0.0128
print('Variacao da razao de carga util total:')
print(DlamT)

ef = 100 * DlamT / 0.0128
print('Variacao percentual na eficiencia:')
print(ef)
