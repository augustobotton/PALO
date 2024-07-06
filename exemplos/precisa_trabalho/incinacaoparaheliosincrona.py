import numpy as np
# Constantes
mu=3.986004418E14 # m^3/s^2

# Dados do problema
J2=0.0010826 # Segunda constante de Jeffery da Terra
Re=6378.14E3 # m - Raio equatorial da Terra
a=6700E3 # m - semi eixo maior da órbita
e=0.01 # excentricidade
# Dado indireto do problema - velocidade hélio síncrona (rotação da Terra ao redor do Sol)
OMpm=2*np.pi/365.25 # rad por dia
OMpm=OMpm/(24*60*60) # rad/s
# Velocidade angular da órbita
n=np.sqrt(mu/(a**3))

print(n)

# Parâmetro da órbita
p=a*(1-e**2)
print(p)

# Inclinação necessária para a órbita hélio síncrona

inc=np.arccos(-OMpm*(2/(3*n*J2))*(p/Re)**2)

print("Inclinação para obter a órbita hélio síncrona do exercício")

print("i = ",inc*180/np.pi,"°")