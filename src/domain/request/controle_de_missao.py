import numpy as np
from scipy.integrate import solve_ivp

from src.domain.aerodinamica.aerodinamica_N_estagios import aerodinamica_multiplos_estagios
from src.domain.modelos.foguete.dinamica_foguete import dinamica_foguete
from src.domain.modelos.foguete.propulsao_N_estagios import propulsao_N_estagios
from src.domain.modelos.planeta.atmosfera.ModeloAtmosferico import ModeloAtmosferico
from src.domain.orbitalUtils.RvelPolar2RvelRet import RvelPolar2RvelRet
from src.domain.orbitalUtils.Vrel2Vine import Vrel2Vine
from src.domain.orbitalUtils.det_orbita import det_orbita
from src.domain.orbitalUtils.long_ECEF2ECI import long_ECEF2ECI
from src.domain.request import parametros

# Inicializa os parâmetros
raio_equatorial_terrestre = parametros.raio_equatorial
velocidade_inercial_de_rotacao_da_terra = parametros.velocidade_inercial_de_rotação_da_terra
mut = parametros.mut
J2 = parametros.J2
J3 = parametros.J3
J4 = parametros.J4
g = parametros.gravidade_padrao_nivel_do_mar
lc = parametros.lc
delta_temperatura_atm = parametros.dT
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
r0 = raio_equatorial_terrestre + h0

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
TF = 25000
# v0 = float(input('Informe o valor inicial da velocidade relativa (m/s): '))
v0 = 1
# phi0 = float(input('Informe a condição inicial do ângulo de elevação (graus): '))
phi0 = 76.8
phi0 = phi0 * np.pi / 180

y = np.cos(ingso) / np.cos(delta0)

if np.abs(y) > 1:
    print('Não eh possivel atingir a inclinacao a partir da latitude inicial. Calculando a menor possivel')
    y = np.sign(y)

Ai_f = np.arcsin(y)
rpgto = raio_equatorial_terrestre + 250e3
agto = (agso + rpgto) / 2
vigto = np.sqrt(mut * (2 / rpgto - 1 / agto))
A0 = np.arctan(
    np.tan(Ai_f) - (rpgto * velocidade_inercial_de_rotacao_da_terra * np.cos(delta0)) / (vigto * np.cos(Ai_f)))

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
    h[i] = X[i, 3] - raio_equatorial_terrestre
    r = X[i, 3]
    delta[i] = X[i, 4]
    lon[i] = X[i, 5]

    # Força propulsiva, massa e ângulos
    ft[i], m[i], mu[i], epsl[i] = propulsao_N_estagios(t[i], X[i])

    # Parâmetros atmosféricos
    modelo_atmosferico = ModeloAtmosferico()
    T[i], _, _, rho[i], _, M[i], _, _, Kn, _, _, R = modelo_atmosferico.calcula(h[i], V[i], lc,
                                                                                delta_temperatura_atm)

    # Forças aerodinâmicas
    D[i], _, _ = aerodinamica_multiplos_estagios(t[i], V[i], h[i], M[i], Kn, T[i], rho[i], R)

    q[i] = 0.5 * rho[i] * V[i] ** 2  # Pressão dinâmica

    # Coordenadas da velocidade inercial no referencial LVLH
    Vi[i], phii[i], Ai[i] = Vrel2Vine(V[i], phi[i], A[i], velocidade_inercial_de_rotacao_da_terra, r, delta[i])

    # Longitude celeste
    longc[i] = long_ECEF2ECI(t[i], lon[i], velocidade_inercial_de_rotacao_da_terra, tg)

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
P = 2 * np.pi * np.sqrt((raio_equatorial_terrestre + hfq[0][0]) ** 3 / mut)  # Período da órbita obtida
print('*** Parametros da Orbita Obtida ***')
print('Velocidade no momento da insercao orbital (km/s):')
print(Vfq[0][0] / 1e3)
print('Altitude no momento da insercao orbital (km):')
print(hfq[0][0] / 1e3)
print('Distancia radial no momento da insercao orbital (km):')
print((hfq[0][0] + raio_equatorial_terrestre) / 1e3)
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
print((rp - raio_equatorial_terrestre) / 1e3)
print('Altitude do apogeu (km):')
print((ra - raio_equatorial_terrestre) / 1e3)
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
    DVgso / (impulso_especico_por_estagio[
                 2] * g))  # Massa de propelente requerida para circularizacao da orbita (kg)
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
