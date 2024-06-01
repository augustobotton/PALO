import numpy as np
import math
from scipy.integrate import solve_ivp
import numpy as np
from aerodinamica_N_estagios import aerodinamica_N_estagios
from propulsao_N_estagios import propulsao_N_estagios
from atm_padrao import atm_padrao
import long_ECEF2ECI
import RvelPolar2RvelRet
import Vrel2Vine
from det_orbita import det_orbita
import matplotlib.pyplot as plt
from grav_axisimetrico import grav_axissimetrico


def dinamica_foguete(t, X):
    
    global  we, Re, lc, dT, h0, l_trilho
    
    # Vetor de estado
    V, A, phi, r, delta = X[0], X[1], X[2], X[3], X[4]
    if V < 0:
        V = 0  # Evita velocidade negativa
    # Funcao para calculo da massa e tracao em funcao do tempo
    ft, m = propulsao_N_estagios(t, ti, tq, ts, Isp, mp, ms, m0, g)
    # Os angulos de apontamento da tubeira sao nulos
    epsl = 0
    mu = 0
    # Funcao para calculo do modelo atmosferico
    h = r - Re  # Altitude
    T, Tm, p, rho, ainf, M, mu, Pr, Kn, d, Re = atm_padrao(h, V, lc)
    # Calculo do modelo aerodinamico
    D, fy, L = aerodinamica_N_estagios(t, V, h, M, Kn, T, rho, Re,ts, Sr, fc)
    # Calculo da gravidade
    gc, gd = grav_axissimetrico(r, delta)
    # Equacoes de cinematica de translacao
    rp = V * np.sin(phi)
    deltap = (V / r) * np.cos(phi) * np.cos(A)
    lonp = (V * np.cos(phi) * np.sin(A)) / (r * np.cos(delta))
    # Equacoes de dinamica de translacao
    Vp = (1 / m) * (ft * np.cos(epsl) * np.cos(mu) - D - m * gc * np.sin(phi) + m * gd * np.cos(phi) * np.cos(A) -
                    m * we**2 * r * np.cos(delta) * (np.cos(phi) * np.cos(A) * np.sin(delta) - np.sin(phi) * np.cos(delta)))
    Ap = (1 / (m * V * np.cos(phi))) * (m * (V**2 / r) * np.cos(phi)**2 * np.sin(A) * np.tan(delta) + ft * np.sin(mu) + fy - m * gd * np.sin(A) +
                                         m * we**2 * r * np.sin(A) * np.sin(delta) * np.cos(delta) - 2 * m * we * V * (np.sin(phi) * np.cos(A) * np.cos(delta) - np.cos(phi) * np.sin(delta)))
    phip = (1 / (m * V)) * (m * (V**2 / r) * np.cos(phi) + ft * np.sin(epsl) * np.cos(mu) + L - m * gc * np.cos(phi) - m * gd * np.sin(phi) * np.cos(A) +
                            m * we**2 * r * np.cos(delta) * (np.sin(phi) * np.cos(A) * np.sin(delta) + np.cos(phi) * np.cos(delta)) + 2 * m * we * V * np.sin(A) * np.cos(delta))
    # Saturacao da altitude
    if h < 0:  # Altitude negativa nao e permitida
        # Mantem as derivadas nulas
        rp = 0
        deltap = 0
        lonp = 0
        Vp = 0
        Ap = 0
        phip = 0
    # Modela o trilho de lancamento
    H = h - h0  # Altura
    if H <= l_trilho and t <= 10:  # Verifica se a altura eh menor que l_trilho nos primeiros segundos da simulacao
        Ap = 0
        phip = 0  # Anula as derivadas dos angulos de orientacao da velocidade
    # Derivada do vetor de estado
    Xp = [Vp, Ap, phip, rp, deltap, lonp]
    return Xp
   
global Re, we, mut, tg, J2, J3, J4, g, lc, dT, Sr, fc, ms, m0, mp, ti, tq, ts, Isp, h0, l_trilho
    
# Mass parameters for structural and payload
ms = np.array([284, 320, 55])
# Adjusting the first two stages
ds = 0.1 * ms[:2]
ms[:2] = ms[:2] - ds
# Structural mass in the third stage
mo = 400
mL = 5
# Propulsive parameters
Isp = np.array([261, 261, 350])
mp = np.array([677+ds[0], 898+ds[1], mo - ms[2] - mL])
ti = np.array([0, 15, 13.5 + 185])
tq = np.array([13.5, 44, ti[2] + 30])
ts = np.array([13.5, 59, tq[1] + 5])

# Aerodynamic and environmental parameters
fc = 1.28
S1 = np.pi*(0.557/2)**2
S2 = np.pi*(0.557/2)**2
S3 = np.pi*(0.46/2)**2
Sp = np.pi*(0.46/2)**2
lt = 12.6
l2 = lt - 3.214
l4 = 0.5
l3 = l2 - 3.294 - l4
f2 = (l2/lt)*0.5 + 0.5
f3 = (l3/lt)*0.5 + 0.5
f4 = (l4/lt)*0.5 + 0.5
Sr = np.array([S1, S2*f2, S3*f3, Sp*f4])
lc, dT = 0.5, 10

# Earth Parameters - axis symmetric model (WGS-84)
Re = 6378.1370e3
we = 7.2921150e-5
g = 9.80665
mut = 3.986004418e14
J2 = 0.00108263
J3 = -0.00000254
J4 = -0.00000161
tg = 0

# Initial Conditions - Alcantara space center (from Google Maps)
h0 = 0
delta0 = -2.3267844*math.pi/180
lon0 = -44.4111042*math.pi/180 - 0.03*math.pi/180
l_trilho = 10

# Calculated parameters
m0 = np.sum(mp) + np.sum(ms) + mL
r0 = Re + h0

# Simplified study by rocket equation
sigma = ms / (ms + mp)
m01 = m0
m02 = m01 - mp[0] - ms[0]
m03 = m02 - mp[1] - ms[1]
lamb = np.array([m02/m01, m03/m02, mL/m03])
lambL = np.prod(lamb)
ve = g * Isp
Dv = -np.sum(ve * np.log(sigma + (1-sigma) * lamb))
Dv3 = -ve[2] * np.log(sigma[2] + (1-sigma[2]) * lamb[2])

# Display data
print(f'Rocket reference area with first stage (m^2): {Sr[0]}')
print(f'Rocket reference area with second stage (m^2): {Sr[1]}')
print(f'Rocket reference area with third stage (m^2): {Sr[2]}')
print(f'Payload reference area (m^2): {Sr[3]}')
print(f'Initial mass before burning each stage - kg: {m0}')
print(f'Structural ratios: {sigma}')
print(f'Payload ratios: {lamb}')
print(f'Total payload ratio: {lambL}')
print(f'Exhaust velocities - m/s: {ve}')
print(f'Total ideal velocity impulse - m/s: {Dv}')
print(f'Ideal velocity impulse impressed by the third stage - m/s: {Dv3}')

# User input data
TF = float(input('Informe o tempo da simulacao (s): '))
v0 = float(input('Informe o valor inicial da velocidade relativa (m/s): '))
phi0 = float(input('Informe a condicao inicial do angulo de elevacao (graus): '))
phi0 = phi0 * math.pi / 180
in_ = float(input('Informe a inclinacao da orbita desejada (graus): '))
in_ = in_ * math.pi / 180

# Perform feasibility test using the orbit inclination and initial latitude
y = math.cos(in_) / math.cos(delta0)
if abs(y) > 1:
    print('Nao eh possivel atingir a inclinacao a partir da latitude inicial. Calculando a menor possivel')
y = math.copysign(1, y)

Ai_f = math.asin(y)  # Final condition of inertial velocity azimuth

# Estimate of the initial condition of relative velocity azimuth
rref = Re + 100e3  # Radius of a circular reference orbit with 100 km altitude
viref = math.sqrt(mut / rref)  # Velocity of a reference orbit with 100 km altitude
A0 = math.atan(math.tan(Ai_f) - (rref * we * math.cos(delta0)) / (viref * math.cos(Ai_f)))

print('Condicao final de azimute de velocidade inercial (°): ', Ai_f * 180 / math.pi)
print('Condicao inicial de azimute de velocidade relativa (°): ', A0 * 180 / math.pi)



# Initial condition
X0 = [v0, A0, phi0, r0, delta0, lon0]
  

t_s = (0,TF)
    
sol = solve_ivp(dinamica_foguete,t_s, y0=X0)

print(sol.success)6


# Post processing
t = sol.t
X = sol.y
N = len(t)  # Number of time instances

 # Initializing arrays for the following variables:
V = np.zeros(N)
A = np.zeros(N)
phi = np.zeros(N)
h = np.zeros(N)
delta = np.zeros(N)
lon = np.zeros(N)
m = np.zeros(N)
ft = np.zeros(N)
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
in_ = np.zeros(N)
om = np.zeros(N)
R0 = np.zeros((N, 3))

for i in range(N):
    V[i] = X[0, i]    # Relative speed
    A[i] = X[1, i]    # Azimuth of relative velocity
    phi[i] = X[2, i]  # Elevation of relative velocity
    h[i] = X[3, i] - Re  # Altitude
    r = X[3, i]  # Radial distance
    delta[i] = X[4, i]  # Latitude
    lon[i] = X[5, i]  # Longitude in the fixed reference to the planet
    ft[i], m[i] = propulsao_N_estagios(t[i])  # Propulsive force and mass
    T[i], _, rho[i], _, M[i], _, _, Kn, _, _ = atm_padrao(h[i], V[i], lc, dT)  # Atmospheric parameters
    D[i], _, _ = aerodinamica_N_estagios(t[i], V[i], h[i], M[i], Kn, T[i], rho[i])  # Aerodynamic forces
    q[i] = 0.5 * rho[i] * V[i]**2  # Dynamic pressure
    # Coordinates of inertial speed in the local horizontal reference
    Vi[i], phii[i], Ai[i] = Vrel2Vine(V[i], phi[i], A[i], we, r, delta[i])
    # Celestial longitude
    longc[i] = long_ECEF2ECI(t[i], lon[i], we, tg)
    # Specific orbital energy
    ee[i] = Vi[i]**2 / 2 - mut / r
    # Position and inertial speed in the celestial ECI reference
    rc0, vc0 = RvelPolar2RvelRet(Vi[i], Ai[i], phii[i], r, delta[i], longc[i])
    R0[i, :] = np.transpose(rc0)
    
    # Orbital parameters
    par_orb = det_orbita(t[i], rc0, vc0, mut)
    a[i], e[i], tau[i], OM[i], in_[i], om[i] = par_orb


# Orbit analysis
  # Generation of vectors for plotting
# Semi-major axis of a reference circular orbit with 100 km altitude
ar = np.ones(N) * rref
# Ve    locity of a reference orbit with 100km altitude
Vir = np.ones(N) * viref
# Specific energy of a reference orbit with a=100km
eer = -mut / (2 * ar)
# Al    titude and inertial speed at the end of third stage burn
for i in range(N):
    if t[i] > tq[2]:
        break
ifq = i - 1
tfq = t[ifq]  # Time of the end of third stage burn
Vfq = np.ones(N) * Vi[ifq]  # Inertial speed at the end of third stage burn
hfq = np.ones(N) * h[ifq]  # Altitude at the end of third stage burn
P = 2 * np.pi * np.sqrt((Re + hfq[0])**3 / mut)  # Orbit period obtained

print('Orbit period obtained (min): ')
print(P / 60)

rp = a[ifq] * (1 - e[ifq])  # Periapsis radius
ra = a[ifq] * (1 + e[ifq])  # Apoapsis radius

print('Periapsis radius (km): ')
print(rp / 1e3)

print('Apoapsis radius (km): ')
print(ra / 1e3)

print('Periapsis altitude (km): ')
print((rp - Re) / 1e3)

print('Apoapsis altitude (km): ')
print((ra - Re) / 1e3)


# Plots
fig1, axs = plt.subplots(2, 3)
axs[0, 0].plot(t, V, linewidth=2)
axs[0, 0].grid(True)
axs[0, 0].set(xlabel='t (s)', ylabel='V (m/s)')

axs[0, 1].plot(t, A * 180 / np.pi, linewidth=2)
axs[0, 1].grid(True)
axs[0, 1].set(xlabel='t (s)', ylabel='A (°)')

axs[0, 2].plot(t, phi * 180 / np.pi, linewidth=2)
axs[0, 2].grid(True)
axs[0, 2].set(xlabel='t (s)', ylabel='φ (°)')

axs[1, 0].plot(t, h / 1e3, linewidth=2)
axs[1, 0].plot(t, hfq / 1e3, '--')
axs[1, 0].plot(tfq, hfq[0] / 1e3, '*')
axs[1, 0].grid(True)
axs[1, 0].set(xlabel='t (s)', ylabel='h (km)')
axs[1, 0].legend(['altitude', 'altitude at end of third stage burn'])

axs[1, 1].plot(t, delta * 180 / np.pi, linewidth=2)
axs[1, 1].grid(True)
axs[1, 1].set(xlabel='t (s)', ylabel='δ (°)')

axs[1, 2].plot(t, lon * 180 / np.pi, linewidth=2)
axs[1, 2].grid(True)
axs[1, 2].set(xlabel='t (s)', ylabel='l(°)')









  

 