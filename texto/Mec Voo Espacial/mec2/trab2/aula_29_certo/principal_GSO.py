import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from dinamica_foguete import dinamica_foguete
from scipy.integrate import solve_ivp
from parametrosglobais import *

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
X0 = [V0, A0, phi0, r0, delta0, lon0]
# Parâmetros para a função de integração
# Simulação
t0 = 0
options = {'rtol': 1e-8, 'atol': 1e-10, 'max_step': 0.5}
resposta_sim = solve_ivp(dinamica_foguete,(t0,Tempo),y0=X0,**options)
t=resposta_sim.t
X = resposta_sim.y