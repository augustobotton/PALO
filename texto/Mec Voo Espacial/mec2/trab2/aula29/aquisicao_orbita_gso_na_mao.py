from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
from solucao_din_foguete import solucao_din_foguete
from atm_padrao import atm_padrao
from Vrel2Vine import Vrel2Vine
from long_ECEF2ECI import long_ECEF2ECI
from RvelPolar2RvelRet import RvelPolar2RvelRet
from scipy.integrate import solve_ivp

# Parâmetros propulsivos
Isp = np.array([251, 271, 315])  # s - Impulso específico dos estagios
mp = np.array([5.5262e+04, 11058, 243.6, 0])  # kg - Massa de propelente dos estagios
Tq1 = 62  # s - NAO MUDA - DADO DOS MOTORES DO 1° ESTAGIO
Tq2 = 64.62  # s - NAO MUDA - DADO DOS MOTORES DO 2° ESTAGIO
Tq3 = 301  # s - TEMPO DE QUEIMA DO 3° ESTAGIO SE ELE IGNITASSE SO UMA VEZ

# Parâmetros de massa estrutural e de carga util
ms = np.array([7750, 1367, 64.7544])  # kg - Massa estrutural dos estagios
mL = 13  # kg - Massa da carga util

# Parâmetros aerodinâmicos e ambientais
fc = 1.28  # Fator de correção do arrasto a partir de dados de tunel de vento
S1 = 4.6*5/3  # m^2 - Area aproximada da seção transversal do primeiro estagio
S2 = 1.5  # m^2 - Area aproximada da seção longitudinal do segundo estagio
S3 = 1.5  # m^2 - Area aproximada da seção longitudinal do terceiro estagio
SL = 1.5  # m^2 -  Area aproximada da seção longitudinal da carga util
lt = 7.33 + 7.1 + 6.28  # m - Comprimento total
l2 = 7.1 + 6.28  #  Comprimento sem o primeiro estagio
l3 = 6.28  #  Comprimento sem o segundo estagio
l4 = 1  # Comprimento da carga util
f2 = (l2/lt)*0.5 + 0.5  # Fator de correção do segundo estagio
f3 = (l3/lt)*0.5 + 0.5  # Fator de correção do terceiro estagio
f4 = (l4/lt)*0.5 + 0.5  # Fator de correção da carga util
Sr = np.array([S1, S2*f2, S3*f3, SL*f4])  # Vetor de areas de referencia para calculo do arrasto
lc = 1.5  # Comprimento característico - diametro dos estagios 2 e superiores
dT = 10  # K - Delta T em relação à atmosfera padrão (que é 15°C no nível do mar)

# Parâmetros da Terra - modelo axis simetrico (WGS-84)
Requat = 6378.1370e3  # m - Raio equatorial da Terra
we = 7.2921150e-5  # (rad/s) - Velocidade inercial de rotação da Terra
g = 9.80665  # m/s^2 - aceleracao da gravidade padrao ao nivel do mar
mut = 3.986004418e14  # m3.s^-2
J2 = 0.00108263  # Constante de Jeffery
J3 = -0.00000254  # Constante de Jeffery
J4 = -0.00000161  # Constante de Jeffery
tg = 0  # s - Tempo em que o meridiano de referência tem longitude celeste nula

# Condições iniciais - Centro espacial de Alcântara (do Google Maps)
h0 = 0  # m - Altitude da base de lançamento
delta0 = -2.3267844 * np.pi / 180  # rad - Latitude inicial
lon0 = -44.4111042 * np.pi / 180  # rad - Longitude inicial
l_trilho = lt  # m - igual ao comprimento total do foguete

# Parâmetros da órbita desejada
ingso = 5 * np.pi / 180  # Inclinação
agso = 42.164140e6  # m
vgso = np.sqrt(mut / agso)

# PARAMETROS PROPULSIVOS E TEMPORAIS A DETERMINAR 
# SAO DEFINIDOS PARA PROPICIAR A INSERCAO ORBITAL
Ts1 = 2  # s
Ts2 = 2  # s
Ts3 = 2  # s
TEq2 = 5  # s
TEq3 = 700 # s
Tq31 = 262
Tq32 = 39
mp31 = mp[2] * Tq31 / Tq3
mp32 = mp[2] * Tq32 / Tq3
mp3 = mp[2]
ti = np.zeros(5)
tq = np.zeros(5)
ts = np.zeros(5)
ti[0] = 0  # s - Tempo da ignicao do primeiro estagio
tq[0] = ti[0] + Tq1  # s - Tempo do fim da queima do est�gio 1
ts[0] = tq[0] + Ts1  # Tempo da separacao do primeiro estagio
ti[1] = ts[0] + TEq2  # s - Tempo da ignicao do segundo estagio
tq[1] = ti[1] + Tq2  # s - Tempo do fim da queima do est�gio 2
ts[1] = tq[1] + Ts2  # Tempo da separacao do segundo estagio
ti[2] = ts[1] + TEq3  # s - Tempo da primeira ignicao do terceiro estagio
tq[2] = ti[2] + Tq31  # s - Tempo do fim da primeira queima do est�gio 3
ti[3] = 1e10  # s - Tempo da segunda ignicao do terceiro estagio
tq[3] = ti[3] + Tq32  # s - Tempo do fim da segunda queima do est�gio 3
ts[2] = tq[3] + Ts3  # Tempo da separacao do terceiro estagio
mp[2] =  mp31 
mp[3] =  mp32 # Complementa o vetor de massas de propelente

# Parametros globais para procurar o apogeu da orbita de transferencia
sinalPhii = 0
achouApogeu = 0

# Parametros calculados
m0 = np.sum(mp) + np.sum(ms) + mL  # Massa inicial do foguete
r0 = Requat + h0  # Distancia radial inicial

# Estudo simplificado pela equacao de foguete
mpx = (mp[0],mp[1], mp3)  # Razão estrutural do primeiro e segundo estágios
print(mpx)
sigma = ms / (ms + mpx)
m01 = m0  # Massa total na decolagem
m02 = ms[1] + mpx[1] + ms[2] + mpx[2] + mL  # Massa total na ignição do segundo estágio
m03 = ms[2] + mpx[2] + mL  # Massa total na ignição do terceiro estágio
lamb = np.zeros(3)
lamb[0] = m02 / m01  # Razão de carga útil do primeiro estágio
lamb[1] = m03 / m02  # Razão de carga útil do segundo estágio
lamb[2] = mL / m03  # Razão de carga útil do terceiro estágio
lambL = np.prod(lamb)  # Razão de carga útil total
ve = g * Isp  # Velocidade de exaustão
Dv = -np.sum(ve * np.log(sigma + (1 - sigma) * lamb))  # Delta v ideal da configuração original

# Mostra dados na tela
print(f'Area de referencia do foguete com primeiro estagio (m^2): {Sr[0]}')
print(f'Area de referencia do foguete com segundo estagio (m^2): {Sr[1]}')
print(f'Area de referencia do foguete com terceiro estagio (m^2): {Sr[2]}')
print(f'Area de referencia da carga util (m^2): {Sr[3]}')
print(f'Massa inicial antes da queima do primeiro estagio - kg: {m01}')
print(f'Massa inicial antes da queima do segundo estagio - kg: {m02}')
print(f'Massa inicial antes da queima do terceiro estagio - kg: {m03}')
print(f'Massa da carga util - kg: {mL}')
print(f'Razoes estruturais: {sigma}')
print(f'Razoes de carga util: {lamb}')
print(f'Velocidades de exaustao - m/s: {ve}')
print(f'Razao de carga util total: {lambL}')
print(f'Impulso de velocidade total ideal - m/s: {Dv}')



Tempo = 5000 #float(input('Informe o tempo da simulacao (s): '))
V0 = 1 #float(input('Informe o valor inicial da velocidade relativa (m/s): '))
phi0 = 84 #float(input('Informe a condicao inicial do angulo de elevacao (graus): '))
phi0 = phi0 * np.pi / 180
y_t = np.cos(ingso) / np.cos(delta0)  # Faz um teste de factibilidade usando a inclinacao da orbita e a latitude inicial
if abs(y_t) > 1:
    print('Nao eh possivel atingir a inclinacao a partir da latitude inicial. Calculando a menor possivel')
    y = np.sign(y_t)
    
Ai_f = np.arcsin(y_t)  # Condicao final do azimute de velocidade inercial

# Estimativa da condicao inicial de azimute de velocidade relativa
rpgto = Requat + 250e3  # Apogeu de uma orbita de transferencia com 250 km de altitude
agto = (agso + rpgto) / 2  # Semi eixo maior de uma orbita de transferencia com 250 km de altitude
vigto = np.sqrt(mut * (2 / rpgto - 1 / agto))  # Velocidade de uma orbita de transferencia com 250km de altitude
A0 = np.arctan(np.tan(Ai_f) - (rpgto * we * np.cos(delta0)) / (vigto * np.cos(Ai_f)))
  
print(f'Condicao final de azimute de velocidade inercial (grau): {Ai_f * 180 / np.pi}')
print(f'Condicao inicial de azimute de velocidade relativa (grau): {A0 * 180 / np.pi}')
    
# Simulação
# Condição inicial
print(V0)
X0 = [V0, A0, phi0, r0, delta0, lon0]
# Parâmetros para a função de integração
# Simulação
print(X0)
t0 = 0
options = {'rtol': 1e-8, 'atol': 1e-10, 'max_step': 0.5}
resposta_sim = solve_ivp(solucao_din_foguete,(t0,Tempo),y0=X0,**options)
t=resposta_sim.t
X = resposta_sim.y


# Pós-processamento
# Parâmetros da órbita GSO requerida
# Velocidade orbital - passada como variável global para a função de dinâmica,
# que deve  calcular o impulso de velocidade de circularização 
print('*** Orbita GSO requerida ***')
print('Raio da orbita GSO (km)', agso/1e3)
print('Velocidade da orbita GSO (km/s)', vgso/1e3)

# Cálculo de outras variáveis
# Número de instantes de tempo
# Magnitude, azimute e elevação da velocidade relativa

N = len(t)
V = resposta_sim.y[0]
A = resposta_sim.y[1]
phi = resposta_sim.y[2]
# Altitude, latitude e longitude no referencial fixo ao planeta

r = resposta_sim.y[3]
h = np.zeros(N)
delta = resposta_sim.y[4]
lon = resposta_sim.y[5]
m = np.zeros(N)  # Massa
ft = np.zeros(N)  # Força propulsiva
mu = np.zeros(N)
epsl = np.zeros(N)  # Ângulos propulsivos
D = np.zeros(N)  # Força de arrasto
q = np.zeros(N)  # Pressão dinâmica
Mach = np.zeros(N)  # Número de Mach
T = np.zeros(N)  # Temperatura
rho = np.zeros(N)  # Densidade
Vi = np.zeros(N)  # Magnitude da velocidade inercial
phii = np.zeros(N)  # Elevação da velocidade inercial
Ai = np.zeros(N)  # Azimute da velocidade inercial
longc = np.zeros(N)  # Longitude celeste
ee = np.zeros(N)  # Energia específica
a = np.zeros(N)  # Semi eixo maior da órbita
e = np.zeros(N)  # Excentricidade da órbita
tau = np.zeros(N)  # Tempo de perigeu
OM = np.zeros(N)  # Ascenção reta do nodo ascendente
in_ = np.zeros(N)  # Inclinação da órbita
om = np.zeros(N)  # Argumento de perigeu
R0 = np.zeros((N, 3))  # Posição no referencial ECI
import propulsao_N_estagios as pn
import aerodinamica_N_estagios as an
for i in range(N):
    
    
    # Posicao no referencial PCPF
    h[i], r[i]- Requat
    # Forca propulsiva, massa e angulos
    ft[i], m[i], mu[i], epsl[i] = pn.propulsao_N_estagios(t[i], V[i], A[i], phi[i], r[i], delta[i],ti, tq, ts, Isp, mp, ms, m0, g, mL, we, Requat )
    # Parametros atmosfericos
    T[i], _, _, rho[i], _, Mach[i], _, _, Kn, _, _, R = atm_padrao(h[i], V[i], lc, dT)
    # Forcas aerodinamicas
    D[i], _, _ = an.aerodinamica_N_estagios(t[i], V[i], h[i], Mach[i], Kn, T[i], rho[i], R,fc,ts,Sr)
    # Pressao dinamica
    q[i] = 0.5 * rho[i] * V[i]**2
    # Coordenadas da velocidade inercial no referencial LVLH
    Vi[i], phii[i], Ai[i] = Vrel2Vine(V[i], phi[i], A[i], we, r[i], delta[i])
    # Longitude celeste
    longc[i] = long_ECEF2ECI(t[i], lon[i], we, tg)
    # Energia especifica da orbita
    ee[i] = Vi[i]**2 / 2 - mut / r[i]
    # Posicao e velocidade inercial no referencial ICP
    rc0, vc0 = RvelPolar2RvelRet(Vi[i], Ai[i], phii[i], r[i], delta[i], longc[i])
    R0[i, :] = rc0.T
    # Elementos orbitais
    par_orb = det_orbita(t[i], rc0, vc0, mut)
    a[i], e[i], tau[i], OM[i], in_[i], om[i] = par_orb[0], par_orb[1], par_orb[2], par_orb[3], par_orb[4], par_orb[5]

# Analise de orbita
# Altitude e velocidade inercial no fim da queima do terceiro estagio
for i in range(N):
    if t[i] > tq[2]:  
        break

ifq = i - 2  # Python's index is 0-based
# Tempo do fim da queima do terceiro estagio
tfq = t[ifq]
# Velocidade inercial no fim da queima do terceiro estagio
Vfq = np.full(N, Vi[ifq])
# Altitude no fim da queima do terceiro estagio
hfq = np.full(N, h[ifq])
# Periodo da orbita obtida
P = 2 * np.pi * np.sqrt((Requat + hfq[0]) ** 3 / mut)
print('*** Parametros da Orbita Obtida ***')
print('Velocidade no momento da insercao orbital (km/s)', Vfq[0] / 1e3)
print('Altitude no momento da insercao orbital (km)', hfq[0] / 1e3)
print('Distancia radial no momento da insercao orbital (km)', (hfq[0] + Requat) / 1e3)
print('Semi eixo maior (km)', a[ifq] / 1e3)
print('Periodo(min): ', P / 60)
# Raio do perigeu
rp = a[ifq] * (1 - e[ifq])
# Raio do apogeu
ra = a[ifq] * (1 + e[ifq])
print('Raio do perigeu (km): ', rp / 1e3)
print('Raio do apogeu (km): ', ra / 1e3)
print('Altitude do perigeu (km): ', (rp - Requat) / 1e3)
print('Altitude do apogeu (km): ', (ra - Requat) / 1e3)

# Orbita de transferencia geossincrona (GTO) desejada
print('*** Parametros da Orbita GTO requerida ***')
print('Perigeu da orbita GTO requerida (km)')
rpgto = rp
print(rpgto / 1e3)
print('Apogeu da orbita GTO requerida (km)')
ragto = agso
print(ragto / 1e3)
print('Semi eixo maior da orbita GTO requerida (km)')
agto = (ragto + rpgto) / 2
print(agto / 1e3)
print('Velocidade de perigeu da orbita GTO requerida (km/s)')
vpgto = np.sqrt(mut * (2 / rpgto - 1 / agto))
print(vpgto / 1e3)
print('Velocidade de apogeu da orbita GTO requerida (km/s)')
vagto = np.sqrt(mut * (2 / ragto - 1 / agto))
print(vagto / 1e3)

# Geracao de vetores para tracar grafico
ar = np.full(N, agto)  # Semi eixo maior da orbita GTO requerida
Vir = np.full(N, vpgto)  # Velocidade de perigeu da orbita GTO requerida
eer = -mut / (2 * ar)  # Energia especifica da orbita GTO requerida
eegso = np.full(N, -mut / (2 * agso))  # Energia especifica da orbita GSO requerida

# Tempos de operacao do propulsor do terceiro estagio
print('Tempo de espera para disparo do propulsor do 3º estagio apos a separacao do 2� (s)', TEq3)
print('Duracao do primeiro disparo do motor do 3º estagio (s)', Tq31)
print('Duracao do segundo disparo do motor do 3º estagio (s)', Tq32)
print('Momento do segundo disparo do motor do 3ºestagio (s)', ti[3])  # Indexing is 0-based in python
print('Impulso de velocidade requerido para circularizacao da orbita (km/s)')
DVgso = vgso - vagto
print(DVgso / 1e3)
print('Massa de propelente requerida para circularizacao da orbita (kg)')
# Massa de propelente necessaria
mp32 = (m[ifq] * np.exp(DVgso / (Isp[2] * g)) - m[ifq]) / np.exp(DVgso / (Isp[2] * g))
print(mp32)
print('Massa de propelente disponivel para o 3� disparo (kg)', mp3 - mp31)
print('****** PARAMETROS DA ORBITA FINAL ******')
print('Periodo (min)')
P = 2 * np.pi * np.sqrt(a[-1] ** 3 / mut)
print(P / 60)
print('Semi eixo maior (km)', a[-1])
print('Excentricidade', e[-1])
print('Inclinacao (�)', in_[-1] * 180 / np.pi)

# Figure 1
plt.figure()
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
plt.ylabel('A ( )')

plt.subplot(233)
plt.plot(t, phi * 180 / np.pi, linewidth=2)
plt.plot(tfq, phi[ifq] * 180 / np.pi, '')
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('\phi ( )')

plt.subplot(234)
plt.plot(t, h / 1e3, linewidth=2)
plt.plot(t, hfq / 1e3, '--')
plt.plot(tfq, hfq[0] / 1e3, '*')
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('h (km)')
plt.legend(['altitude', 'altitude no fim da queima do 3  estagio'])

plt.subplot(235)
plt.plot(t, delta * 180 / np.pi, linewidth=2)
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('\delta ( )')

plt.subplot(236)
plt.plot(t, lon * 180 / np.pi, linewidth=2)
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('l( )')

# Figure 2
plt.figure()
plt.subplot(221)
plt.plot(t, Vi, linewidth=2)
plt.plot(t, Vir, '--')
plt.plot(t, Vfq, '-.')
plt.plot(tfq, Vfq[0], '*')
plt.grid(True)
plt.xlabel('t (s)')
plt.ylabel('V_i (m/s)')
plt.legend(['Velocidade inercial', 'Velocidade de perigeu da orbita GTO requerida', 'Velocidade no fim da queima do terceiro estagio'])

plt.subplot(222)
plt.plot(t, Ai * 180 / np.pi, linewidth=2)
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('A_i ( )')

plt.subplot(223)
plt.plot(t, phii * 180 / np.pi, linewidth=2)
plt.plot(tfq, phii[ifq] * 180 / np.pi, '')
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('\phi_i ( )')

plt.subplot(224)
plt.plot(t, longc * 180 / np.pi, linewidth=2)
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('\lambda ( )')

# Figure 3
plt.figure()
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
plt.ylabel('\mu ( )')

plt.subplot(224)
plt.plot(t, epsl * 180 / np.pi, linewidth=2)
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('\epsilon ( )')

# Figure 4
plt.figure()
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
plt.plot(t, Mach, linewidth=2)
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('M (-)')

plt.subplot(325)
plt.plot(t, T - 273.15, linewidth=2)
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('T ( C)')

plt.subplot(326)
plt.plot(t, rho, linewidth=2)
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('\rho (kg/m^3)')

# Figure 5
plt.figure()
plt.subplot(311)
plt.plot(t, ee, linewidth=2)
plt.plot(t, eer, '--', linewidth=2)
plt.plot(t, eegso, '--', linewidth=2)
plt.grid(True)
plt.xlabel('t (s)')
plt.ylabel('\epsilon (J/kg)')
plt.legend(['Energia especifica', 'Energia especifica da orbita GTO requerida', 'Energia especifica da orbita GSO requerida'])

plt.subplot(334)
plt.plot(t, a / 1e3, linewidth=2)
plt.plot(t, ar / 1e3, '--')
plt.plot(t, np.full((N,), Requat) / 1e3, '-.')
plt.grid(True)
plt.xlabel('t (s)')
plt.ylabel('a (km)')
plt.legend(['Semi eixo maior', 'Semi eixo maior da orbita GTO requerida', 'Raio da Terra'])

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
plt.ylabel('\tau (s)')

plt.subplot(337)
plt.plot(t, OM * 180 / np.pi, linewidth=2)
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('\Omega ( )')

plt.subplot(338)
plt.plot(t, in_ * 180 / np.pi, linewidth=2)  # renamed 'in' to 'in_' as 'in' is a reserved word in Python
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('i ( )')

plt.subplot(339)
plt.plot(t, om * 180 / np.pi, linewidth=2)
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('\omega (grau)')


def earth_sphere(num_points=100):
    # Definindo as coordenadas para o globo terrestre
    u = np.linspace(0, 2 * np.pi, num_points)
    v = np.linspace(0, np.pi, num_points)
    x = 6371 * np.outer(np.cos(u), np.sin(v))
    y = 6371 * np.outer(np.sin(u), np.sin(v))
    z = 6371 * np.outer(np.ones(np.size(u)), np.cos(v))

    return x, y, z


def plot_trajectory_on_earth(latitudes, longitudes, altitudes):
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')

    x, y, z = earth_sphere()

    # Desenho do globo terrestre
    ax.plot_surface(x, y, z, color='dodgerblue', alpha=0.3, linewidth=0)

    # Conversão das coordenadas do trajeto para coordenadas cartesianas
    xs = (6371 + altitudes) * np.cos(latitudes) * np.cos(longitudes)
    ys = (6371 + altitudes) * np.cos(latitudes) * np.sin(longitudes)
    zs = (6371 + altitudes) * np.sin(latitudes)

    # Plotando a trajetória
    ax.plot(xs, ys, zs, color='darkorange')

    ax.set_xlabel('X [km]')
    ax.set_ylabel('Y [km]')
    ax.set_zlabel('Z [km]')

    plt.show()

plot_trajectory_on_earth(np.rad2deg(delta),np.rad2deg(lon),h)



plt.show()