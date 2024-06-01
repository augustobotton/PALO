import numpy as np
import matplotlib.pyplot as plt
from solucao_din_foguete import solucao_din_foguete
from aerodinamica_N_estagios import aerodinamica_N_estagios
from atm_padrao import atm_padrao
from propulsao_N_estagios import propulsao_N_estagios
from Vrel2Vine import Vrel2Vine
from det_orbita import det_orbita
from long_ECEF2ECI import long_ECEF2ECI
from RvelPolar2RvelRet import RvelPolar2RvelRet
# Script para simular um voo de sondagem de veiculo de dois estagios com
# carga util.
# Eh utilizado o foguete de insercao orbital geosincrono conceitual
# (FIOGSC)
# de referencias indicadas em aula.

## Parametros propulsivos
Isp = np.array([251, 271, 315])     # s - Specific impulse of the first stage (5xS50), second stage (1xS50), and second stage (RD843)
mp = np.array([5.5262e+04, 11058, 243.6])     # kg - Propellant mass of the first stage (5xS50), second stage (1xS50), and second stage (RD843)
Tq1, Tq2, Tq3 = 62, 64.62, 301
ti = np.zeros(3) # s - Ignition time of the first stage
tq = np.zeros(3)    # s - End time of burning stage 1
ts = np.zeros(3) # Time of separation of the first stage

# Time calculations
tq[0] = ti[0] + Tq1
ts[0] = tq[0] + 1
ti[1] = ts[0] + 1
tq[1] = ti[1] + Tq2
ts[1] = tq[1] + 1
ti[2] = ts[1] + 1
tq[2] = ti[2] + Tq3
ts[2] = tq[2] + 1

# Structural and payload mass parameters
ms = np.array([7750, 1367, 64.7544]) # kg - Structural mass of the first stage, second stage, and second stage
mL = 13 # kg - Payload mass

## Parametros aerodinamicos e ambientais
# Fator de correcao do arrasto a partir de dados de tunel de vento
fc=1.28
# As areas de referencia sao das secoes transversais
S1=4.6*5/3 # m^2 - Area aproximada da secao transversal do primeiro estagio
S2=1.5 # m^2 - Area aproximada da secao longitudinal do segundo estagio
S3=1.5 # m^2 - Area aproximada da secao longitudinal do terceiro estagio
SL=1.5 # m^2 -  Area aproximada da secao longitudinal da carga util
# Correcoes para levar em conta a area molhada, feitas em funcao do
# comprimento. Assume-se que a correcao soh ocorre em 50% da area de
# referencia, supondo uma influencia de 50% da area molhada no arrasto
# total
lt=7.33+7.1+6.28 # m - Comprimento total
l2=7.1+6.28    #  Comprimento sem o primeiro estagio
l3=6.28    #  Comprimento sem o segundo estagio
l4=1 # Comprimento da carga util
f2=(l2/lt)*0.5+0.5 # Fator de correcao do segundo estagio
f3=(l3/lt)*0.5+0.5 # Fator de correcao do terceiro estagio
f4=(l4/lt)*0.5+0.5 # Fator de correcao da carga util
# Vetor de areas de referencia para calculo do arrasto
Sr=[S1,
   S2*f2,
   S3*f3,
   SL*f4]
lc=1.5 # Comprimento caracteristico - diametro dos estagios 2 e superiores
dT=10  # K - Delta T em relação à atmosfera padr o (que   15 C no n vel do mar)

# Parametros da Terra - modelo axis simetrico (WGS-84)
Requat = 6378.1370e3  # m - Raio equatorial da Terra
we = 7.2921150e-5  # (rad/s) - Velocidade inercial de rotação da Terra
g = 9.80665  # m/s^2 - aceleracao da gravidade padrao ao nivel do mar
mut = 3.986004418e14  # m3.s^-2
# Constantes de Jeffery
J2 = 0.00108263
J3 = -0.00000254
J4 = -0.00000161
tg = 0   # s - Tempo em que o meridiano de referência tem longitude celeste nula

# Condicoes iniciais - Centro espacial de Alcantara (do Google Maps)
h0 = 0  # m - Altitude da base de lancamento
delta0 = -2.3267844 * np.pi/180  # rad - Latitude inicial
lon0 = -44.4111042 * np.pi/180  # rad - Longitude inicial
# Comprimento do trilho de lancamento
l_trilho = lt  # m - igual ao comprimento total do foguete

# Dados fornecidos no teclado pelo usuario
Tempo = 600 #float(input('Informe o tempo da simulacao (s): '))
V0 = 5  #float(input('Informe o valor inicial da velocidade relativa (m/s): '))
A0 = 85 #float(input('Informe o azimute inicial da velocidade relativa (graus): '))
A0 = A0 * np.pi/180  # Azimute da velocidade relativa em t=0
phi0 =85 # float(input('Informe a elevacao inicial da velocidade relativa (graus): '))
phi0 = phi0 * np.pi/180  # Elevação da velocidade relativa em t=0

# Parametros calculados
# Massa inicial do foguete
m0 = sum(mp) + sum(ms) + mL
# Distancia radial inicial
r0 = Requat + h0

# Estudo simplificado pela equacao de foguete
# Razão estrutural do primeiro e segundo estágios
sigma = np.divide(ms, (ms + mp))
# Massa total na decolagem
m01 = m0
# Massa total na ignição do segundo estágio
m02 = ms[1] + mp[1] + ms[2] + mp[2] + mL
# Massa total na ignição do terceiro estágio
m03 = ms[2] + mp[2] + mL
# Razão de carga útil do primeiro, segundo e terceiro estágios
lamb = [m02/m01, m03/m02, mL/m03]
# Razão de carga útil total
lambL = np.prod(lamb)

# Velocidade de exaustão
ve = g * np.array(Isp)
# Delta v ideal da configuração original
Dv = -sum(ve * np.log(np.add(sigma, np.multiply((1 - sigma), lamb))))

# Mostra dados na tela
print('Area de referencia do foguete com primeiro estagio (m^2):', Sr[0])
print('Area de referencia do foguete com segundo estagio (m^2):', Sr[1])
print('Area de referencia do foguete com terceiro estagio (m^2):', Sr[2])
print('Area de referencia da carga util (m^2):', Sr[3])
print('Massa inicial antes da queima do primeiro estagio - kg', m01)
print('Massa inicial antes da queima do segundo estagio - kg', m02)
print('Massa inicial antes da queima do terceiro estagio - kg', m03)
print('Massa da carga util - kg', mL)
print('Razoes estruturais', sigma)
print('Razoes de carga util', lamb)
print('Velocidades de exaustao - m/s', ve)
print('Razao de carga útil total', lambL)
print('Impulso de velocidade total ideal - m/s', Dv)

# Simulação
# Condição inicial
X0 = [V0, A0, phi0, r0, delta0, lon0]
# Parâmetros para a função de integração
# Simulação

t0 = 0
resposta_sim =solucao_din_foguete(t0,Tempo, X0, we, Requat, lc, dT, h0, l_trilho, ti, ts, tq, Isp, mp, ms, m0, g,fc,Sr)
t=resposta_sim.t

N = len(t)
V = resposta_sim.y[0]
A = resposta_sim.y[1]
phi = resposta_sim.y[2]
r = resposta_sim.y[3]
h = np.zeros(N)
delta = resposta_sim.y[4]
lon = resposta_sim.y[5]
m = np.zeros(N)
ft = np.zeros(N)
D = np.zeros(N)
q = np.zeros(N)
Mach = np.zeros(N)
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
in_ = np.zeros(N)
om = np.zeros(N)
R0 = np.zeros((N, 3))


for i in range(N):
    # Magnitude, azimute e elevação da velocidade relativa
    
    # Posição na referencial do PCPF
    h[i] = r[i] - Requat
            
    # Força propulsiva e massa    
    ft[i], m[i] = propulsao_N_estagios(t[i],ti, ts, tq, Isp, mp, ms, m0, g)  
    
    # Atmospheric parameters
    T[i], _, _, rho[i], _, Mach[i], _, _, Kn, _, _, R = atm_padrao(h[i], V[i], lc, dT)
    
    # Aerodynamic forces
    D[i], _, _ = aerodinamica_N_estagios(t[i], V[i], h[i], Mach[i], Kn, T[i], rho[i], R,fc,ts,Sr)   
    
    # Dynamic pressure
    q[i] = 0.5 * rho[i] * V[i]**2 
    
    # Inertial speed coordinates in LVLH reference
    
    Vi[i], phii[i], Ai[i] = Vrel2Vine(V[i], phi[i], A[i], we, r[i], delta[i])
    
    # Celestial longitude
    longc[i] = long_ECEF2ECI(t[i], lon[i], we, tg)
    
    # Specific orbital energy
    ee[i] = Vi[i]**2 / 2 - mut / r[i]

    # Position and inertial velocity in ICP reference
    rc0, vc0 = RvelPolar2RvelRet(Vi[i], Ai[i], phii[i], r[i], delta[i], longc[i])
    R0[i, :] = np.transpose(rc0)

    # Orbital elements
    par_orb = det_orbita(t[i], rc0, vc0, mut)
    a[i], e[i], tau[i] = par_orb[0], par_orb[1], par_orb[2]
    OM[i], in_[i], om[i] = par_orb[3], par_orb[4], par_orb[5]




# Orbit analysis
n_gso = we  # Mean angular speed of geosynchronous orbit
a_gso = (mut/n_gso**2)**(1/3)  # Semi major axis of geosynchronous orbit
viref = np.sqrt(mut/a_gso)  # Speed of a circular GSO reference orbit

# Generation of vectors to plot graph
ar = a_gso * np.ones(N)  # Semi major axis of a GSO reference circular orbit
Vir = viref * np.ones(N)  # Speed of a GSO reference orbit
eer = -mut / (2 * ar)  # Specific energy of a GSO reference orbit

# Altitude and inertial speed at the end of the third stage burn
for i in range(N):
    if t[i] > tq[2]:
        break

ifq = i - 1
tfq = t[ifq]  # Time of the end of the third stage burn
Vfq = Vi[ifq] * np.ones(N)  # Inertial speed at the end of the third stage burn
hfq = h[ifq] * np.ones(N)  # Altitude at the end of the third stage burn
P = 2 * np.pi * np.sqrt(((Requat + hfq[0]) ** 3 )/ mut)  # Period of the obtained orbit

print('Period of the obtained orbit (min): ')
print(P / 60)

rp = a[ifq] * (1 - e[ifq])  # Periapsis radius
ra = a[ifq] * (1 + e[ifq])  # Apoapsis radius

print('Periapsis radius (km): ')
print(rp / 1e3)
print('Apoapsis radius (km): ')
print(ra / 1e3)
print('Periapsis altitude (km): ')
print((rp - Requat) / 1e3)
print('Apoapsis altitude (km): ')
print((ra - Requat) / 1e3)


fig1, axs = plt.subplots(2, 3, figsize=(15,10))

axs[0, 0].plot(t, V, linewidth=2)
axs[0, 0].set(xlabel='t (s)', ylabel='V (m/s)')
axs[0, 0].grid(True)

axs[0, 1].plot(t, A * 180 / np.pi, linewidth=2)
axs[0, 1].set(xlabel='t (s)', ylabel='A (º)')
axs[0, 1].grid(True)

axs[0, 2].plot(t, phi * 180 / np.pi, linewidth=2)
axs[0, 2].set(xlabel='t (s)', ylabel='\phi (º)')
axs[0, 2].grid(True)

axs[1, 0].plot(t, h / 1e3, linewidth=2)
axs[1, 0].plot(t, hfq / 1e3, '--')
axs[1, 0].plot(tfq, hfq[0] / 1e3, '*')
axs[1, 0].set(xlabel='t (s)', ylabel='h (km)')
axs[1, 0].grid(True)
axs[1, 0].legend(['altitude', 'altitude at the end of the 3rd stage burn'])

axs[1, 1].plot(t, delta * 180 / np.pi, linewidth=2)
axs[1, 1].set(xlabel='t (s)', ylabel='\delta (º)')
axs[1, 1].grid(True)

axs[1, 2].plot(t, lon * 180 / np.pi, linewidth=2)
axs[1, 2].set(xlabel='t (s)', ylabel='l(º)')
axs[1, 2].grid(True)

fig2, axs2 = plt.subplots(2, 2, figsize=(15,10))

axs2[0, 0].plot(t, Vi, linewidth=2)
axs2[0, 0].plot(t, Vir, '--')
axs2[0, 0].plot(t, Vfq, '-.')
axs2[0, 0].plot(tfq, Vfq[0], '*')
axs2[0, 0].set(xlabel='t (s)', ylabel='V_i (m/s)')
axs2[0, 0].legend(['Velocidade inercial', 'Velocidade de orbita circular de referencia GSO',
            'Velocidade no fim da queima do terceiro estagio'])
axs2[0, 0].grid(True)

axs2[0, 1].plot(t, Ai * 180 / np.pi, linewidth=2)
axs2[0, 1].set(xlabel='t (s)', ylabel='A_i (º)')
axs2[0, 1].grid(True)

axs2[1, 0].plot(t, phii * 180 / np.pi, linewidth=2)
axs2[1, 0].set(xlabel='t (s)', ylabel='\phi_i (º)')
axs2[1, 0].grid(True)

axs2[1, 1].plot(t, longc * 180 / np.pi, linewidth=2)
axs2[1, 1].set(xlabel='t (s)', ylabel='\lambda (º)')
axs2[1, 1].grid(True)

fig3, axs3 = plt.subplots(2, 1, figsize=(10,10))

axs3[0].plot(t, ft, linewidth=2)
axs3[0].set(xlabel='t (s)', ylabel='f_t (N)')
axs3[0].grid(True)

axs3[1].plot(t, m, linewidth=2)
axs3[1].set(xlabel='t (s)', ylabel='m (kg)')
axs3[1].grid(True)


fig4, axs4 = plt.subplots(3, 2, figsize=(15,15))

axs4[0, 0].plot(t, D, linewidth=2)
axs4[0, 0].set(xlabel='t (s)', ylabel='D (N)')
axs4[0, 0].grid(True)

axs4[0, 0].plot(t, V, linewidth=2)
axs4[0, 0].set(xlabel='t (s)', ylabel='D (N)')
axs4[0, 0].grid(True)

axs4[1, 0].plot(t, q, linewidth=2)
axs4[1, 0].set(xlabel='t (s)', ylabel='q (N/m^2)')
axs4[1, 0].grid(True)

axs4[1, 1].plot(t, Mach, linewidth=2)
axs4[1, 1].set(xlabel='t (s)', ylabel='M (-)')
axs4[1, 1].grid(True)

axs4[2, 0].plot(t, T-273.15, linewidth=2)
axs4[2, 0].set(xlabel='t (s)', ylabel='T (ºC)')
axs4[2, 0].grid(True)

axs4[2, 1].plot(t, rho, linewidth=2)
axs4[2, 1].set(xlabel='t (s)', ylabel='\rho (kg/m^3)')
axs4[2, 1].grid(True)


fig5, axs5 = plt.subplots(3, 3, figsize=(15,15))

axs5[0, 0].plot(t, ee, t, eer, '--', linewidth=2)
axs5[0, 0].set(xlabel='t (s)', ylabel='\epsilon (J/kg)')
axs5[0, 0].legend(['Energia especifica','Energia especifica de orbita GSO'])
axs5[0, 0].grid(True)

axs5[1, 0].plot(t, a/1e3, linewidth=2)
axs5[1, 0].plot(t, ar/1e3, '--')
axs5[1, 0].plot(t, np.ones(N)*Requat/1e3, '-.')
axs5[1, 0].set(xlabel='t (s)', ylabel='a (km)')
axs5[1, 0].legend(['Semi eixo maior', 'Semi eixo maior de orbita GSO', 'Raio da Terra'])
axs5[1, 0].grid(True)

axs5[1, 1].plot(t, e, linewidth=2)
axs5[1, 1].set(xlabel='t (s)', ylabel='e (-)')
axs5[1, 1].grid(True)

axs5[1, 2].plot(t, tau, linewidth=2)
axs5[1, 2].set(xlabel='t (s)', ylabel='\tau (s)')
axs5[1, 2].grid(True)

axs5[2, 0].plot(t, OM*180/np.pi, linewidth=2)
axs5[2, 0].set(xlabel='t (s)', ylabel='\Omega (º)')
axs5[2, 0].grid(True)

axs5[2, 1].plot(t, in_*180/np.pi, linewidth=2)
axs5[2, 1].set(xlabel='t (s)', ylabel='i (º)')
axs5[2, 1].grid(True)

axs5[2, 2].plot(t, om*180/np.pi, linewidth=2)
axs5[2, 2].set(xlabel='t (s)', ylabel='\omega (º)')
axs5[2, 2].grid(True)

plt.tight_layout()
plt.show()