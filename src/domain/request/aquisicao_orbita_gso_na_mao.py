import os

import numpy as np
from src.domain.modelos.planeta import ModeloAtmosferico
from scipy.integrate import solve_ivp
from src.domain.modelos.foguete import ModeloAerodinamico

import parametros
from src.domain.modelos.foguete.dinamica_foguete import dinamica_foguete
from src.domain.modelos.foguete import ModeloPropulsivo
from src.domain.orbitalUtils.Converte import *
from src.domain.orbitalUtils.det_orbita import det_orbita


global Re, we, mut, J2, J3, J4, g, lc, dT, Sr, fator_correcao, massa_carga_util
global ms, m0, mp, ti, tq, ts, Isp, h0, l_trilho, tg, agso, Tq3, Tq31, Tq32, Ts3, vgso, mp3


def clear_console():
    os.system('cls' if os.name == 'nt' else 'clear')


# Parâmetros propulsivos
Isp = np.array([251, 271, 315])  # s - Impulso específico dos estágios
parametros.impulso_especifico_por_estagio = Isp
mp = np.array([5.5262e4, 11058, 0, 0])  # kg - Massa de propelente dos estágios
parametros.mp = mp
mp3 = 243.6  # kg - Massa de propelente do terceiro estágio (RD843)
Tq1 = 62  # s - DADO DOS MOTORES DO 1° ESTÁGIO
Tq2 = 64.62  # s - DADO DOS MOTORES DO 2° ESTÁGIO
Tq3 = 301  # s - TEMPO DE QUEIMA DO 3° ESTÁGIO SE ELE IGNITASSE SÓ UMA VEZ

# Parâmetros de massa estrutural e de carga útil
ms = np.array([7750, 1367, 64.7544])  # kg - Massa estrutural dos estágios
parametros.massa_estrutural_por_estagio = ms
massa_carga_util = 13  # kg - Massa da carga útil
parametros.massa_de_carga_util = massa_carga_util

# Parâmetros aerodinâmicos e ambientais
fator_correcao = 1.28  # Fator de correção do arrasto
parametros.fator_correcao_arrasto = fator_correcao
area_secao_transversal_primeiro_estagio = 4.6 * 5 / 3
area_secao_transversal_segundo_estagio = 1.5  # m^2 - Area aproximada da secao transversal do segundo estagio
area_secao_transversal_terceiro_estagio = 1.5  # m^2 - Area aproximada da secao transversal do terceiro estagio
SL = 1.5  # m^2 - Area aproximada da secao transversal da carga útil

lt = 7.33 + 7.1 + 6.28  # m - Comprimento total
l2 = 7.1 + 6.28  # Comprimento sem o primeiro estagio
l3 = 6.28  # Comprimento sem o segundo estagio

l4 = 1  # Comprimento da carga util
f2 = (l2 / lt) * 0.5 + 0.5  # Fator de correcao do segundo estagio
f3 = (l3 / lt) * 0.5 + 0.5  # Fator de correcao do terceiro estagio
f4 = (l4 / lt) * 0.5 + 0.5  # Fator de correcao da carga util

Sr = [area_secao_transversal_primeiro_estagio, area_secao_transversal_segundo_estagio * f2,
      area_secao_transversal_terceiro_estagio * f3, SL * f4]
parametros.area_de_referencia = Sr

lc = 1.5  # Comprimento característico - diâmetro dos estágios 2 e superiores
parametros.lc = lc
dT = 10  # K - Delta T em relação à atmosfera padrão
parametros.dT = dT
Re = 6378.1370e3  # m - Raio equatorial da Terra
parametros.raio_equatorial = Re
we = 7.2921150e-5  # (rad/s) - Velocidade inercial de rotação da Terra
parametros.velocidade_inercial_de_rotação_da_terra = we
g = 9.80665  # m/s^2 - aceleração da gravidade ao nível do mar
parametros.gravidade_padrao_nivel_do_mar = g
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

ingso = 5 * np.pi / 180  # Inclinacao
agso = 42.164140e6  # m
vgso = np.sqrt(mut / agso)

# Parametros propulsivos e temporais a determinar
Ts1 = 2
Ts2 = 2
Ts3 = 2
TEq2 = 5  # s
TEq3 = 1360  # s

Tq31 = 262
Tq32 = 39

mp31 = mp3 * Tq31 / Tq3
mp32 = mp3 * Tq32 / Tq3

ti = np.array([0, 0, 0, 0])
tq = np.array([0, 0, 0, 0])
ts = np.array([0, 0, 0])

ti[0] = 0
tq[0] = ti[0] + Tq1
ts[0] = tq[0] + Ts1
ti[1] = ts[0] + TEq2
tq[1] = ti[1] + Tq2
ts[1] = tq[1] + Ts2
ti[2] = ts[1] + TEq3
tq[2] = ti[2] + Tq31

ti[3] = 1e6
tq[3] = ti[3] + Tq32
ts[2] = tq[3] + Ts2

mp[2] = mp31
mp[3] = mp32

parametros.ti = ti
parametros.tq = tq
parametros.tempo_limite_separacao = ts

sinalPhii = 0
parametros.sinalPhii = sinalPhii
achouApogeu = 0
parametros.achouApogeu = achouApogeu

m0 = np.sum(mp) + np.sum(ms) + massa_carga_util
parametros.m0 = m0
r0 = Re + h0

# Estudo simplificado pela equação de foguete
mpx = [mp[0], mp[1], mp3]
ms_mpx_sum = ms + mpx
sigma = ms / ms_mpx_sum

m01 = m0
m02 = ms[1] + mpx[1] + ms[2] + mpx[2] + massa_carga_util
m03 = ms[2] + mpx[2] + massa_carga_util

lamb0 = m02 / m01
lamb1 = m03 / m02
lamb2 = massa_carga_util / m03
lamb = np.array([lamb0, lamb1, lamb2])

lambL = np.prod(lamb)

ve = g * Isp

Dv = -np.sum(ve * np.log(sigma + (1 - sigma) * lamb))

# Mostra dados na tela
print('Area de referencia do foguete com primeiro estagio (m^2):', Sr[0])
print('Area de referencia do foguete com segundo estagio (m^2):', Sr[1])
print('Area de referencia do foguete com terceiro estagio (m^2):', Sr[2])
print('Area de referencia da carga util (m^2):', Sr[3])
print('Massa inicial antes da queima do primeiro estagio - kg:', m01)
print('Massa inicial antes da queima do segundo estagio - kg:', m02)
print('Massa inicial antes da queima do terceiro estagio - kg:', m03)
print('Massa da carga util - kg:', massa_carga_util)
print('Razoes estruturais:', sigma)
print('Razoes de carga util:', lamb)
print('Velocidades de exaustao - m/s:', ve)
print('Razao de carga útil total:', lambL)
print('Impulso de velocidade total ideal - m/s:', Dv)

# Simulação
simula = 1

while simula == 1:
    TF = 6000  # float(input('Informe o tempo da simulação (s): '))
    v0 = 2  # float(input('Informe o valor inicial da velocidade relativa (m/s): '))
    phi0 = 87  # float(input('Informe a condição inicial do ângulo de elevação (graus): '))
    phi0 = phi0 * np.pi / 180

    y = np.cos(ingso) / np.cos(delta0)
    print(y)

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

    X0 = np.array([v0, A0, phi0, r0, delta0, lon0])
    options = {
        'rtol': 1e-8,
        'atol': 1e-10,
        'max_step': 0.5
    }

    # Solve the differential equations
    sol = solve_ivp(dinamica_foguete, [0, TF], X0, method='LSODA', **options)

    t = sol.t
    X = sol.y
    X = np.transpose(X)
    break

# Parâmetros da órbita GSO requerida
print('*** Órbita GSO requerida ***')
print('Raio da órbita GSO (km):', agso / 1e3)
print('Velocidade da órbita GSO (km/s):', vgso / 1e3)

# Cálculo de outras variáveis
N = len(t)  # Número de instantes de tempo

# Inicialização dos vetores
V = np.zeros(N)
A = np.zeros(N)
phi = np.zeros(N)
h = np.zeros(N)
delta = np.zeros(N)
lon = np.zeros(N)
m = np.zeros(N)
ft = np.zeros(N)
mu = np.zeros(N)
epsl = np.zeros(N)
D = np.zeros(N)
q = np.zeros(N)
M = np.zeros(N)
T = np.zeros(N)
rho = np.zeros(N)
Vi = np.zeros(N)
phii = np.zeros(N)
Ai = np.zeros(N)
longc = np.zeros(N)
ee = np.zeros(N)
a = np.zeros(N)
e = np.zeros(N)
tau = np.zeros(N)
OM = np.zeros(N)
inclinacao = np.zeros(N)
om = np.zeros(N)
R0 = np.zeros((N, 3))

for i in range(0, N):
    # Magnitude, azimute e elevação da velocidade relativa
    V[i] = X[i, 0]
    A[i] = X[i, 1]
    phi[i] = X[i, 2] - 1e-4

    # Posição no referencial PCPF
    h[i] = X[i, 3] - Re
    r = X[i, 3]
    delta[i] = X[i, 4]
    lon[i] = X[i, 5]

    # Força propulsiva, massa e ângulos
    ft[i], m[i], mu[i], epsl[i] = propulsao_N_estagios(t[i], X[i, :])

    # Parâmetros atmosféricos
    T[i], _, _, rho[i], _, M[i], _, _, Kn, _, _, R = atm_padrao(h[i], V[i], lc, dT)

    # Forças aerodinâmicas
    D[i], _, _ = aerodinamica_multiplos_estagios(t[i], V[i], h[i], M[i], Kn, T[i], rho[i], R)

    q[i] = 0.5 * rho[i] * V[i] ** 2  # Pressão dinâmica

    # Coordenadas da velocidade inercial no referencial LVLH
    Vi[i], phii[i], Ai[i] = Vrel2Vine(V[i], phi[i], A[i], we, r, delta[i])

    # Longitude celeste
    longc[i] = converte_longitudeFixaPlaneta_para_longitude_celeste(t[i], lon[i], we, tg)

    # Energia específica da órbita
    ee[i] = Vi[i] ** 2 / 2 - mut / r

    # Posição e velocidade inercial no referencial ICP
    rc0, vc0 = RvelPolar2RvelRet(Vi[i], Ai[i], phii[i], r, delta[i], longc[i])
    R0[i, :] = rc0

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
ifq = i - 2
tfq = t[ifq]  # Tempo do fim da queima do terceiro estágio
Vfq = Vi[ifq] * np.ones(N)  # Velocidade inercial no fim da queima do terceiro estágio
hfq = h[ifq] * np.ones(N)  # Altitude no fim da queima do terceiro estágio
P = 2 * np.pi * np.sqrt((Re + hfq[0]) ** 3 / mut)  # Período da órbita obtida
print('*** Parametros da Orbita Obtida ***')
print('Velocidade no momento da insercao orbital (km/s):')
print(Vfq[0] / 1e3)
print('Altitude no momento da insercao orbital (km):')
print(hfq[0] / 1e3)
print('Distancia radial no momento da insercao orbital (km):')
print((hfq[0] + Re) / 1e3)
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
ar = agto * np.ones(N)  # Semi eixo maior da orbita GTO requerida
Vir = vpgto * np.ones(N)  # Velocidade de perigeu da orbita GTO requerida
eer = -mut / (2 * ar)  # Energia especifica da orbita GTO requerida
eegso = -mut / (2 * agso) * np.ones(N)  # Energia especifica da orbita GSO requerida
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
mp32 = (m[ifq] * np.exp(DVgso / (Isp[2] * g)) - m[ifq]) / np.exp(
    DVgso / (Isp[2] * g))  # Massa de propelente requerida para circularizacao da orbita (kg)
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
# plt.figure(1)
#
# plt.subplot(231)
# plt.plot(t, V, linewidth=2)
# plt.grid(True)
# plt.axis('tight')
# plt.xlabel('t (s)')
# plt.ylabel('V (m/s)')
#
# plt.subplot(232)
# plt.plot(t, A * 180 / np.pi, linewidth=2)
# plt.grid(True)
# plt.axis('tight')
# plt.xlabel('t (s)')
# plt.ylabel('A (º)')
#
# plt.subplot(233)
# plt.plot(t, phi * 180 / np.pi, linewidth=2)
# plt.plot(tfq, phi[ifq] * 180 / np.pi, '*')
# plt.grid(True)
# plt.axis('tight')
# plt.xlabel('t (s)')
# plt.ylabel('phi (º)')
#
# plt.subplot(234)
# plt.plot(t, h / 1e3, linewidth=2)
# plt.plot(t, hfq / 1e3, '--')
# plt.plot(tfq, hfq[0] / 1e3, '*')
# plt.grid(True)
# plt.axis('tight')
# plt.xlabel('t (s)')
# plt.ylabel('h (km)')
# plt.legend(['altitude', 'altitude no fim da queima do 3º estágio'])
#
# plt.subplot(235)
# plt.plot(t, delta * 180 / np.pi, linewidth=2)
# plt.grid(True)
# plt.axis('tight')
# plt.xlabel('t (s)')
# plt.ylabel('delta (º)')
#
# plt.subplot(236)
# plt.plot(t, lon * 180 / np.pi, linewidth=2)
# plt.grid(True)
# plt.axis('tight')
# plt.xlabel('t (s)')
# plt.ylabel('l(º)')
#
# # Figure 2
# plt.figure(2)
#
# plt.subplot(221)
# plt.plot(t, Vi, linewidth=2)
# plt.plot(t, Vir, '--')
# plt.plot(t, Vfq, '-.')
# plt.plot(tfq, Vfq[0], '*')
# plt.grid(True)
# plt.xlabel('t (s)')
# plt.ylabel('V_i (m/s)')
# plt.legend(['Velocidade inercial', 'Velocidade de perigeu da órbita GTO requerida',
#             'Velocidade no fim da queima do terceiro estágio'])
#
# plt.subplot(222)
# plt.plot(t, Ai * 180 / np.pi, linewidth=2)
# plt.grid(True)
# plt.axis('tight')
# plt.xlabel('t (s)')
# plt.ylabel('A_i (º)')
#
# plt.subplot(223)
# plt.plot(t, phii * 180 / np.pi, linewidth=2)
# plt.plot(tfq, phii[ifq] * 180 / np.pi, '*')
# plt.grid(True)
# plt.axis('tight')
# plt.xlabel('t (s)')
# plt.ylabel('phi_i (º)')
#
# plt.subplot(224)
# plt.plot(t, longc * 180 / np.pi, linewidth=2)
# plt.grid(True)
# plt.axis('tight')
# plt.xlabel('t (s)')
# plt.ylabel('lambda (º)')
#
# # Figure 3
# plt.figure(3)
#
# plt.subplot(221)
# plt.plot(t, ft, linewidth=2)
# plt.grid(True)
# plt.axis('tight')
# plt.xlabel('t (s)')
# plt.ylabel('f_t (N)')
#
# plt.subplot(222)
# plt.plot(t, m, linewidth=2)
# plt.grid(True)
# plt.axis('tight')
# plt.xlabel('t (s)')
# plt.ylabel('m (kg)')
#
# plt.subplot(223)
# plt.plot(t, mu * 180 / np.pi, linewidth=2)
# plt.grid(True)
# plt.axis('tight')
# plt.xlabel('t (s)')
# plt.ylabel('\u03BC (º)')
#
# plt.subplot(224)
# plt.plot(t, epsl * 180 / np.pi, linewidth=2)
# plt.grid(True)
# plt.axis('tight')
# plt.xlabel('t (s)')
# plt.ylabel('\u03B5 (º)')
#
# # Figure 4
# plt.figure(4)
#
# plt.subplot(311)
# plt.plot(t, D, linewidth=2)
# plt.grid(True)
# plt.axis('tight')
# plt.xlabel('t (s)')
# plt.ylabel('D (N)')
#
# plt.subplot(323)
# plt.plot(t, q, linewidth=2)
# plt.grid(True)
# plt.axis('tight')
# plt.xlabel('t (s)')
# plt.ylabel('q (N/m^2)')
#
# plt.subplot(324)
# plt.plot(t, M, linewidth=2)
# plt.grid(True)
# plt.axis('tight')
# plt.xlabel('t (s)')
# plt.ylabel('M (-)')
#
# plt.subplot(325)
# plt.plot(t, T - 273.15, linewidth=2)
# plt.grid(True)
# plt.axis('tight')
# plt.xlabel('t (s)')
# plt.ylabel('T (ºC)')
#
# plt.subplot(326)
# plt.plot(t, rho, linewidth=2)
# plt.grid(True)
# plt.axis('tight')
# plt.xlabel('t (s)')
# plt.ylabel('rho (kg/m^3)')
#
# # Figure 5
# plt.figure(5)
#
# plt.subplot(311)
# plt.plot(t, ee, linewidth=2)
# plt.plot(t, eer, '--', linewidth=2)
# plt.plot(t, eegso, '--', linewidth=2)
# plt.grid(True)
# plt.xlabel('t (s)')
# plt.ylabel('\u03B5 (J/kg)')
# plt.legend(['Energia específica', 'Energia específica da órbita GTO requerida',
#             'Energia específica da órbita GSO requerida'])
#
# plt.subplot(334)
# plt.plot(t, a / 1e3, linewidth=2)
# plt.plot(t, ar / 1e3, '--')
# plt.plot(t, Re * np.ones(N) / 1e3, '-.')
# plt.grid(True)
# plt.xlabel('t (s)')
# plt.ylabel('a (km)')
# plt.legend(['Semi eixo maior', 'Semi eixo maior da órbita GTO requerida', 'Raio da Terra'])
#
# plt.subplot(335)
# plt.plot(t, e, linewidth=2)
# plt.grid(True)
# plt.axis('tight')
# plt.xlabel('t (s)')
# plt.ylabel('e (-)')
#
# plt.subplot(336)
# plt.plot(t, tau, linewidth=2)
# plt.grid(True)
# plt.axis('tight')
# plt.xlabel('t (s)')
# plt.ylabel('\u03C4 (s)')
#
# plt.subplot(337)
# plt.plot(t, OM * 180 / np.pi, linewidth=2)
# plt.grid(True)
# plt.axis('tight')
# plt.xlabel('t (s)')
# plt.ylabel('\u03A9 (º)')
#
# plt.subplot(338)
# plt.plot(t, inclinacao * 180 / np.pi, linewidth=2)
# plt.grid(True)
# plt.axis('tight')
# plt.xlabel('t (s)')
# plt.ylabel('i (º)')
#
# plt.subplot(339)
# plt.plot(t, om * 180 / np.pi, linewidth=2)
# plt.grid(True)
# plt.axis('tight')
# plt.xlabel('t (s)')
# plt.ylabel('\u03C9 (º)')
#
# # Figure 6
# fig6 = plt.figure(6)
# ax6 = fig6.add_subplot(111, projection='3d')
# ax6.plot(delta * 180 / np.pi, lon * 180 / np.pi, h, linewidth=2)
# ax6.set_xlabel('\u03B4 (º)')
# ax6.set_ylabel('long (º)')
# ax6.set_zlabel('h')
# ax6.set_xlim(-180, 180)
# ax6.set_ylim(-180, 180)
# ax6.set_zlim(0, np.max(h))
#
# # Figure 7
# fig7 = plt.figure(7)
# ax7 = fig7.add_subplot(111, projection='3d')
# ax7.plot(R0[:, 0] / 1e3, R0[:, 1] / 1e3, R0[:, 2] / 1e3, linewidth=2)
# ax7.set_xlabel('X (km)')
# ax7.set_ylabel('Y (km)')
# ax7.set_zlabel('Z (km)')
# ax7.set_xlim(-np.max(R0[:, 0]) / 1e3, np.max(R0[:, 0]) / 1e3)
# ax7.set_ylim(-np.max(R0[:, 1]) / 1e3, np.max(R0[:, 1]) / 1e3)
# ax7.set_zlim(-np.max(R0[:, 2]) / 1e3, np.max(R0[:, 2]) / 1e3)
# ax7.grid(True)
# ax7.set_aspect('equal')
# ax7.plot_surface(np.zeros((2, 2)), np.zeros((2, 2)), np.zeros((2, 2)), alpha=0)  # Plot empty surface for aspect ratio
#
# # Mostra os gráficos
# plt.show()
#
simula = int(input('Deseja simular novamente? (1) sim (0) nao: '))

# fig = plt.figure(figsize=(12, 8))
# ax = fig.add_subplot(111)
# desenha_mapa_trajetoria([delta0 * 180 / np.pi, 2 * lon0 * 180 / np.pi, h0], traj, ax=ax)
# plt.show()
