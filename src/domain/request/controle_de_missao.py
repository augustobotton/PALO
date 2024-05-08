import numpy as np
from scipy.integrate import solve_ivp
from src.domain.modelos.foguete.aerodinamica_N_estagios import aerodinamica_multiplos_estagios
import parametros
from src.domain.modelos.foguete.dinamica_foguete import dinamica_foguete
from src.domain.modelos.foguete.propulsao_N_estagios import propulsao_N_estagios
from src.domain.modelos.planeta.ModeloAtmosferico import ModeloAtmosferico
from src.domain.orbitalUtils.Converte import *
from src.domain.orbitalUtils.det_orbita import det_orbita


ingso = 2.3267844 * np.pi / 180  # Inclinacao
agso = 42.164140e6  # m
parametros.agso = agso
vgso = np.sqrt(mut / agso)
parametros.vgso = vgso

tempo_de_simulacao = 25000
valor_inicial_velocidade_relativa = 1
angulo_de_elevacao_inicial = np.deg2rad(78)

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

X0 = [valor_inicial_velocidade_relativa, A0, angulo_de_elevacao_inicial, r0, delta0, lon0]
# X0 = X0.T
options = {
    'rtol': 1e-12,
    'atol': 1e-14,
    # 'max_step': 0.1
}

# Solve the differential equations
# sol = solve_ivp(dinamica_foguete, [0, TF], X0, method='LSODA', **options)
sol = solve_ivp(dinamica_foguete, [0, tempo_de_simulacao], X0, **options)
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
    longc[i] = converte_longitudeFixaPlaneta_para_longitude_celeste(t[i], lon[i], velocidade_inercial_de_rotacao_da_terra, tg)

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
