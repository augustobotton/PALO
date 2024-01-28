import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

import parametros
from RvelPolar2RvelRet import RvelPolar2RvelRet
from Vrel2Vine import Vrel2Vine
from aerodinamica_N_estagios import aerodinamica_multiplos_estagios
from atm_padrao import atm_padrao
from dinamica_foguete import dinamica_foguete
from domain.OrbitalUtils.det_orbita import det_orbita
from long_ECEF2ECI import long_ECEF2ECI
from propulsao_N_estagios import propulsao_N_estagios

# Inicializa os parâmetros
Re = parametros.Re
we = parametros.we
mut = parametros.mut
J2 = parametros.J2
J3 = parametros.J3
J4 = parametros.J4
g = parametros.g
lc = parametros.lc
dT = parametros.dT
Sr = parametros.area_de_referencia
fator_correcao = parametros.fator_correcao_arrasto
massa_de_carga_util = parametros.massa_de_carga_util
massa_estrutural_por_estagio = parametros.massa_estrutural_por_estagio
m0 = parametros.m0
mp = parametros.mp
ti = parametros.ti
tq = parametros.tq
ts = parametros.tempo_limite_separacao
impulso_especico_por_estagio = parametros.impulso_especifico_por_estagio
h0 = parametros.h0
l_trilho = parametros.comprimento_do_trilho
tg = parametros.tg
agso = parametros.agso
Tq3 = parametros.Tq3
Tq31 = parametros.Tq31
Tq32 = parametros.Tq32
Ts3 = parametros.Ts3
vgso = parametros.vgso
mp3 = parametros.mp3

# Parâmetros propulsivos
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
dT = 10  # K - Delta T em relação à atmosfera padrão
parametros.dT = dT
Re = 6378.1370e3  # m - Raio equatorial da Terra
parametros.Re = Re
we = 7.2921150e-5  # (rad/s) - Velocidade inercial de rotação da Terra
parametros.we = we
g = 9.80665  # m/s^2 - aceleração da gravidade ao nível do mar
parametros.g = g
mut = 3.986004418e14  # m^3/s^-2
parametros.mut = mut
J2 = 0.00108263  # Constante de Jeffery J2
parametros.J2 = J2
J3 = -0.00000254  # Constante de Jeffery J3
parametros.J3 = J3
J4 = -0.00000161  # Constante de Jeffery J4
parametros.J4 = J4
tg = 0  # s - Tempo em que o meridiano de referência tem longitude celeste nula
parametros.tg = tg
# Condições iniciais
h0 = 0  # m - Altitude da base de lancamento
parametros.h0 = h0
delta0 = -2.3267844 * np.pi / 180  # rad - Latitude inicial
lon0 = -44.4111042 * np.pi / 180  # rad - Longitude inicial
l_trilho = lt  # m - igual ao comprimento total do foguete
parametros.comprimento_do_trilho = l_trilho

ingso = 2.3267844 * np.pi / 180  # Inclinacao
agso = 42.164140e6  # m
parametros.agso = agso
vgso = np.sqrt(mut / agso)
parametros.vgso = vgso

# Parametros propulsivos e temporais a determinar
Ts1 = 2
Ts2 = 2
Ts3 = 2
parametros.Ts1 = Ts1
parametros.Ts2 = Ts2
parametros.Ts3 = Ts3
TEq2 = 5  # s
TEq3 = 300  # s

Tq31 = 222.8  # 222.8
parametros.Tq31 = Tq31
Tq32 = Tq3 - Tq31
parametros.Tq32 = Tq32

mp31 = parametros.mp31
mp32 = parametros.mp32

mp31 = mp3 * Tq31 / Tq3
parametros.mp31 = mp31
mp32 = mp3 * Tq32 / Tq3
parametros.mp32 = mp32

ti[0] = 0
parametros.ti[0] = ti[0]
tq[0] = ti[0] + Tq1
parametros.tq[0] = tq[0]
ts[0] = tq[0] + Ts1
parametros.tempo_limite_separacao[0] = ts[0]
ti[1] = ts[0] + TEq2
parametros.ti[1] = ti[1]
tq[1] = ti[1] + Tq2
parametros.tq[1] = tq[1]
ts[1] = tq[1] + Ts2
parametros.tempo_limite_separacao[1] = ts[1]
ti[2] = ts[1] + TEq3
parametros.ti[2] = ti[2]
tq[2] = ti[2] + Tq31
parametros.tq[2] = tq[2]
ti[3] = 1e10
parametros.ti[3] = ti[3]
tq[3] = ti[3] + Tq32
parametros.tq[3] = tq[3]
ts[2] = tq[3] + Ts3
parametros.tempo_limite_separacao[2] = ts[2]

parametros.ti = ti
parametros.tq = tq
parametros.tempo_limite_separacao = ts

mp[2] = mp31
parametros.mp[2] = mp[2]
mp[3] = mp32
parametros.mp[3] = mp[3]
parametros.mp = mp

sinalPhii = parametros.sinalPhii
achouApogeu = parametros.achouApogeu

m0 = np.sum(mp) + np.sum(massa_estrutural_por_estagio) + massa_de_carga_util
parametros.m0 = m0
r0 = Re + h0

# Estudo simplificado pela equação de foguete
mpx = np.array([mp[0], mp[1], mp3])
ms_mpx_sum = massa_estrutural_por_estagio + mpx
sigma = massa_estrutural_por_estagio / ms_mpx_sum

m01 = m0
m02 = massa_estrutural_por_estagio[1] + mpx[1] + massa_estrutural_por_estagio[2] + mpx[2] + massa_de_carga_util
m03 = massa_estrutural_por_estagio[2] + mpx[2] + massa_de_carga_util

lamb0 = m02 / m01
lamb1 = m03 / m02
lamb2 = massa_de_carga_util / m03
lamb = np.array([lamb0, lamb1, lamb2])

lambL = np.prod(lamb)

ve = g * impulso_especico_por_estagio

Dv = -np.sum(ve * np.log(sigma + (1 - sigma) * lamb))

# Mostra dados na tela
print('Area de referencia do foguete com primeiro estagio (m^2):', Sr[0])
print('Area de referencia do foguete com segundo estagio (m^2):', Sr[1])
print('Area de referencia do foguete com terceiro estagio (m^2):', Sr[2])
print('Area de referencia da carga util (m^2):', Sr[3])
print('Massa inicial antes da queima do primeiro estagio - kg:', m01)
print('Massa inicial antes da queima do segundo estagio - kg:', m02)
print('Massa inicial antes da queima do terceiro estagio - kg:', m03)
print('Massa da carga util - kg:', massa_de_carga_util)
print('Razoes estruturais:', sigma)
print('Razoes de carga util:', lamb)
print('Velocidades de exaustao - m/s:', ve)
print('Razao de carga útil total:', lambL)
print('Impulso de velocidade total ideal - m/s:', Dv)

# TF = float(input('Informe o tempo da simulação (s): '))
# v0 = float(input('Informe o valor inicial da velocidade relativa (m/s): '))
# phi0 = float(input('Informe a condição inicial do ângulo de elevação (graus): '))
TF = 25000
v0 = 1
phi0 = 76.8
phi0 = phi0 * np.pi / 180

y = np.cos(ingso) / np.cos(delta0)

if np.abs(y) > 1:
    print('Não eh possivel atingir a inclinacao a partir da latitude inicial. Calculando a menor possivel')
    y = np.sign(y)

Ai_f = np.arcsin(y)
rpgto = Re + 250e3
agto = (agso + rpgto) / 2
vigto = np.sqrt(mut * (2 / rpgto - 1 / agto))
A0 = np.arctan(np.tan(Ai_f) - (rpgto * we * np.cos(delta0)) / (vigto * np.cos(Ai_f)))

print('Condicao final de azimute de velocidade inercial (º):')
print(Ai_f * 180 / np.pi)
print('Condicao inicial de azimute de velocidade relativa (º):')
print(A0 * 180 / np.pi)

X0 = [v0, A0, phi0, r0, delta0, lon0]
# X0 = X0.T
options = {
    'rtol': 1e-12,
    'atol': 1e-14,
    # 'max_step': 0.1
}

# Solve the differential equations
# sol = solve_ivp(dinamica_foguete, [0, TF], X0, method='LSODA', **options)
sol = solve_ivp(dinamica_foguete, [0, TF], X0, **options)
t = sol.t  # Redimensiona t para ser uma matriz de coluna
X = sol.y.T

# t = sol.t
# X = sol.y
# X = np.transpose(X)

# Parâmetros da órbita GSO requerida
print('*** Órbita GSO requerida ***')
print('Raio da órbita GSO (km):', agso / 1e3)
print('Velocidade da órbita GSO (km/s):', vgso / 1e3)

# Cálculo de outras variáveis
N = len(t)  # Número de instantes de tempo

# Inicialização dos vetores
V = np.zeros([N, 1])
A = np.zeros([N, 1])
phi = np.zeros([N, 1])
h = np.zeros([N, 1])
delta = np.zeros([N, 1])
lon = np.zeros([N, 1])
m = np.zeros([N, 1])
ft = np.zeros([N, 1])
mu = np.zeros([N, 1])
epsl = np.zeros([N, 1])
D = np.zeros([N, 1])
q = np.zeros([N, 1])
M = np.zeros([N, 1])
T = np.zeros([N, 1])
rho = np.zeros([N, 1])
Vi = np.zeros([N, 1])
phii = np.zeros(N)
Ai = np.zeros([N, 1])
longc = np.zeros([N, 1])
ee = np.zeros([N, 1])
a = np.zeros([N, 1])
e = np.zeros([N, 1])
tau = np.zeros([N, 1])
OM = np.zeros([N, 1])
inclinacao = np.zeros([N, 1])
om = np.zeros([N, 1])
R0 = np.zeros([N, 3])

for i in range(N):
    # Magnitude, azimute e elevação da velocidade relativa
    V[i] = X[i, 0]
    A[i] = X[i, 1]
    phi[i] = X[i, 2]

    # Posição no referencial PCPF
    h[i] = X[i, 3] - Re
    r = X[i, 3]
    delta[i] = X[i, 4]
    lon[i] = X[i, 5]

    # Força propulsiva, massa e ângulos
    ft[i], m[i], mu[i], epsl[i] = propulsao_N_estagios(t[i], X[i])

    # Parâmetros atmosféricos
    T[i], _, _, rho[i], _, M[i], _, _, Kn, _, _, R = atm_padrao(h[i], V[i], lc, dT)

    # Forças aerodinâmicas
    D[i], _, _ = aerodinamica_multiplos_estagios(t[i], V[i], h[i], M[i], Kn, T[i], rho[i], R)

    q[i] = 0.5 * rho[i] * V[i] ** 2  # Pressão dinâmica

    # Coordenadas da velocidade inercial no referencial LVLH
    Vi[i], phii[i], Ai[i] = Vrel2Vine(V[i], phi[i], A[i], we, r, delta[i])

    # Longitude celeste
    longc[i] = long_ECEF2ECI(t[i], lon[i], we, tg)

    # Energia específica da órbita
    ee[i] = Vi[i] ** 2 / 2 - mut / r

    # Posição e velocidade inercial no referencial ICP
    rc0, vc0 = RvelPolar2RvelRet(Vi[i], Ai[i], phii[i], r, delta[i], longc[i])
    R0[i, :] = rc0.T

    # Elementos orbitais
    par_orb = det_orbita(t[i], rc0, vc0, mut)
    a[i] = par_orb[0]
    e[i] = par_orb[1]
    tau[i] = par_orb[2]
    OM[i] = par_orb[3]
    inclinacao[i] = par_orb[4]
    om[i] = par_orb[5]

# Análise de órbita
for i in range(N):
    if t[i] > tq[2]:
        break

if inclinacao[-1] > 90:
    inclinacao[-1] = 180 - inclinacao[-1]

ifq = i - 2
tfq = t[ifq]  # Tempo do fim da queima do terceiro estágio
Vfq = Vi[ifq] * np.ones([1, N])  # Velocidade inercial no fim da queima do terceiro estágio
hfq = h[ifq] * np.ones([1, N])  # Altitude no fim da queima do terceiro estágio
P = 2 * np.pi * np.sqrt((Re + hfq[0][0]) ** 3 / mut)  # Período da órbita obtida
print('*** Parametros da Orbita Obtida ***')
print('Velocidade no momento da insercao orbital (km/s):')
print(Vfq[0][0] / 1e3)
print('Altitude no momento da insercao orbital (km):')
print(hfq[0][0] / 1e3)
print('Distancia radial no momento da insercao orbital (km):')
print((hfq[0][0] + Re) / 1e3)
print('Semi eixo maior (km):')
print(a[ifq - 1] / 1e3)
print('Periodo (min): ')
print(P / 60)
rp = a[ifq - 1] * (1 - e[ifq - 1])  # Raio do perigeu
ra = a[ifq - 1] * (1 + e[ifq - 1])  # Raio do apogeu
print('Raio do perigeu (km):')
print(rp / 1e3)
print('Raio do apogeu (km):')
print(ra / 1e3)
print('Altitude do perigeu (km):')
print((rp - Re) / 1e3)
print('Altitude do apogeu (km):')
print((ra - Re) / 1e3)
print('*** Parametros da Orbita GTO requerida ***')
print('Perigeu da orbita GTO requerida (km):')
rpgto = rp
print(rpgto / 1e3)
print('Apogeu da orbita GTO requerida (km):')
ragto = agso
print(ragto / 1e3)
print('Semi eixo maior da orbita GTO requerida (km):')
agto = (ragto + rpgto) / 2
print(agto / 1e3)
print('Velocidade de perigeu da orbita GTO requerida (km/s):')
vpgto = np.sqrt(mut * (2 / rpgto - 1 / agto))
print(vpgto / 1e3)
print('Velocidade de apogeu da orbita GTO requerida (km/s):')
vagto = np.sqrt(mut * (2 / ragto - 1 / agto))
print(vagto / 1e3)
ar = agto * np.ones([1, N])  # Semi eixo maior da orbita GTO requerida
Vir = vpgto * np.ones([1, N])  # Velocidade de perigeu da orbita GTO requerida
eer = -mut / (2 * ar)  # Energia especifica da orbita GTO requerida
print('Energia Específica da Órbita GTO Requerida')
# print(eer)
eegso = -mut / (2 * agso) * np.ones([1, N])  # Energia especifica da orbita GSO requerida
print('Tempo de espera para disparo do propulsor do 3º estagio apos a separacao do 2º (s):')
print(TEq3)
print('Duracao do primeiro disparo do motor do 3º estagio (s):')
print(Tq31)
print('Duracao do segundo disparo do motor do 3º estagio (s):')
print(Tq32)
print('Momento do segundo disparo do motor do 3º estagio (s):')
print(ti[3])
DVgso = vgso - vagto  # Impulso de velocidade requerido para circularizacao da orbita (km/s)
print('Impulso de velocidade requerido para circularizacao da orbita (km/s):')
print(DVgso / 1e3)
mp32 = (m[ifq - 1] * np.exp(DVgso / (impulso_especico_por_estagio[2] * g)) - m[ifq - 1]) / np.exp(
    DVgso / (impulso_especico_por_estagio[2] * g))  # Massa de propelente requerida para circularizacao da orbita (kg)
print('Massa de propelente requerida para circularizacao da orbita (kg):')
print(mp32)
print('Massa de propelente disponivel para o 3º disparo (kg):')
print(mp3 - mp31)
print('****** PARAMETROS DA ORBITA FINAL ******')
print('Periodo (min):')
P = 2 * np.pi * np.sqrt(a[-1] ** 3 / mut)
print(P / 60)
print('Semi eixo maior (km):')
print(a[-1])
print('Excentricidade:')
print(e[-1])
print('Inclinacao (º):')
print(inclinacao[-1] * 180 / np.pi)

# Gráficos
# Figure 1
plt.close('all')
plt.figure(1)

plt.subplot(231)
plt.plot(t, V, linewidth=2)
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('V (m/s)')

plt.subplot(232)
plt.plot(t, A * 180 / np.pi, linewidth=2)
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('A (º)')

plt.subplot(233)
plt.plot(t, phi * 180 / np.pi, linewidth=2)
plt.plot(tfq, phi[ifq - 1] * 180 / np.pi, '*')
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('phi (º)')

plt.subplot(234)
plt.plot(t, h / 1e3, linewidth=2)
plt.plot(t, hfq.T / 1e3, '--')
plt.plot(tfq, hfq[0][0] / 1e3, '*')
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('h (km)')
plt.legend(['altitude', 'altitude no fim da queima do 3º estágio'])

plt.subplot(235)
plt.plot(t, delta * 180 / np.pi, linewidth=2)
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('delta (º)')

plt.subplot(236)
plt.plot(t, lon * 180 / np.pi, linewidth=2)
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('l(º)')

# Figure 2
plt.figure(2)

plt.subplot(221)
plt.plot(t, Vi, linewidth=2)
plt.plot(t, Vir.T, '--')
plt.plot(t, Vfq.T, '-.')
plt.plot(tfq, Vfq[0][0], '*')
plt.grid(True)
plt.xlabel('t (s)')
plt.ylabel('V_i (m/s)')
plt.legend(['Velocidade inercial', 'Velocidade de perigeu da órbita GTO requerida',
            'Velocidade no fim da queima do terceiro estágio'])

plt.subplot(222)
plt.plot(t, Ai * 180 / np.pi, linewidth=2)
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('A_i (º)')

plt.subplot(223)
plt.plot(t, phii * 180 / np.pi, linewidth=2)
plt.plot(tfq, phii[ifq - 1] * 180 / np.pi, '*')
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('phi_i (º)')

plt.subplot(224)
plt.plot(t, longc * 180 / np.pi, linewidth=2)
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('lambda (º)')

# Figure 3
plt.figure(3)

plt.subplot(221)
plt.plot(t, ft, linewidth=2)
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('f_t (N)')

plt.subplot(222)
plt.plot(t, m, linewidth=2)
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('m (kg)')

plt.subplot(223)
plt.plot(t, mu * 180 / np.pi, linewidth=2)
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('\u03BC (º)')

plt.subplot(224)
plt.plot(t, epsl * 180 / np.pi, linewidth=2)
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('\u03B5 (º)')

# Figure 4
plt.figure(4)

plt.subplot(311)
plt.plot(t, D, linewidth=2)
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('D (N)')

plt.subplot(323)
plt.plot(t, q, linewidth=2)
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('q (N/m^2)')

plt.subplot(324)
plt.plot(t, M, linewidth=2)
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('M (-)')

plt.subplot(325)
plt.plot(t, T - 273.15, linewidth=2)
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('T (ºC)')

plt.subplot(326)
plt.plot(t, rho, linewidth=2)
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('rho (kg/m^3)')

# Figure 5
plt.figure(5)

plt.subplot(311)
plt.plot(t, ee, linewidth=2)
plt.plot(t, eer.T, '--', linewidth=2)
plt.plot(t, eegso.T, '--', linewidth=2)
plt.grid(True)
plt.xlabel('t (s)')
plt.ylabel('\u03B5 (J/kg)')
plt.legend(['Energia específica', 'Energia específica da órbita GTO requerida',
            'Energia específica da órbita GSO requerida'])

plt.subplot(334)
plt.plot(t, a / 1e3, linewidth=2)
plt.plot(t, ar.T / 1e3, '--')
plt.plot(t, Re * np.ones([N, 1]) / 1e3, '-.')
plt.grid(True)
plt.xlabel('t (s)')
plt.ylabel('a (km)')
plt.legend(['Semi eixo maior', 'Semi eixo maior da órbita GTO requerida', 'Raio da Terra'])

plt.subplot(335)
plt.plot(t, e, linewidth=2)
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('e (-)')

plt.subplot(336)
plt.plot(t, tau, linewidth=2)
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('\u03C4 (s)')

plt.subplot(337)
plt.plot(t, OM * 180 / np.pi, linewidth=2)
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('\u03A9 (º)')

plt.subplot(338)
plt.plot(t, inclinacao * 180 / np.pi, linewidth=2)
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('i (º)')

plt.subplot(339)
plt.plot(t, om * 180 / np.pi, linewidth=2)
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('\u03C9 (º)')

# Figure 6
plt.figure(5)

traj = np.column_stack((delta, lon)) * 180 / np.pi
##desenha_mapa_trajetoria([delta0 * 180 / np.pi, lon0 * 180 / np.pi, h0], traj)
plt.show()

# Figure 7
fig7 = plt.figure(7)
ax = plt.axes(projection="3d")

u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
r = Re / 1e3

x = r * np.outer(np.cos(u), np.sin(v))
y = r * np.outer(np.sin(u), np.sin(v))
z = r * np.outer(np.ones(np.size(u)), np.cos(v))

ax.plot_surface(x, y, z, rstride=4, cstride=4)
ax.plot3D(R0[:, 0] / 1e3, R0[:, 1] / 1e3, R0[:, 2] / 1e3, 'red')
ax = plt.gca()
ax.set_aspect('equal', adjustable='box')
ax.set_xlabel('x (km)')
ax.set_ylabel('y (km)')
ax.set_zlabel('z (km)')

# Mostra os gráficos
plt.show()
