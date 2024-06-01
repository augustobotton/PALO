g = 9.81
# Propulsao
mp1 = 33157 * 5 / 3
ms1 = 4650 * 5 / 3
mp2 = 11058
ms2 = 1367
mp3 = 0.4 * 609
ms3 = (0.21 / (1 - 0.21)) * mp3
Isp1 = 251
Isp2 = 271
Isp3 = 315
F1 = 1317e3 * 5 / 3
F2 = 455e3
F3 = 2.5e3
mp1p = F1 / (g * Isp1)
mp2p = F2 / (g * Isp2)
mp3p = F3 / (g * Isp3)
Tq1 = mp1 / mp1p
Tq2 = mp2 / mp2p
Tq3 = mp3 / mp3p

# Aerodinamica
S1 = 4.6 * 5 / 3
S2 = 1.5
S3 = 1.5
l1 = 7.33
l2 = 7.1
l3 = 6.28

print('Modelo do foguete conceitual')
print('Massas estruturais de cada estagio - kg')
print(ms1)
print(ms2)
print(ms3)
print('Massas de propelente de cada estagio - kg')
print(mp1)
print(mp2)
print(mp3)
print('Impulsos especificos de cada estagio - s')
print(Isp1)
print(Isp2)
print(Isp3)
print('Tempos de queima de cada estagio - s')
print(Tq1)
print(Tq2)
print(Tq3)
