# Importa os pacotes salvos
import perturbacoes
from src.domain.utilidades_mecanica_orbital.orbitalUtils import calculos_orbitais
from src.domain.utilidades_mecanica_orbital.propagacao import propagacao
# Importações de pacotes do Python
import numpy as np
from scipy.integrate import solve_ivp # Solução numérica de equações diferenciais
import matplotlib.pyplot as plt # Gráficos
# Variáveis globais
mu=1.32712440018e20 # m^3/s^2 Parâmetro gravitacional do Sol (primário)
mu3=3.986e14; # m^3/s^2 - Parâmetro gravitacional da Terra (terceiro corpo)
# Vetor de elementos orbitais da Terra com respeito ao sol
at=149597870E3 # m
et= 0.01667
taut=-100*24*60*60 # s
orbt=np.array([at,et,taut,0,0,0])
# Vetor de estado inicial do VE - posição e velocidade
R0=np.array([-27,147.5,0.1])*1E9 # m
V0=np.array([-33,-10,1])*1E3 # m/s
# Função para calcular a dinâmica da massa de prova m2 perturbada pelo terceiro corpo m3
def dinPert3corpo(t,X):
    """
    Entradas:
    t: [s] - tempo
    X: vetor de estado, 6x1
    Saída
    xp: derivada de x
    """
    # Posição e velocidade
    R =np.array(X[0:3]);V=np.array(X[3:6])
    #Distância
    r=np.linalg.norm(R)
    # Posição do terceiro corpo com respeito ao primário
    _,R31,_=propagacao.propagaEliptica(t,mu,orbt)
    # Aceleração perturbativa
    ad= perturbacoes.perturb3Corpo(R, R31, mu3)
    # Cinemática
    Rp=V
    # Dinâmica
    Vp=-mu*R/r**3+ad
    #Saída
    Xp=np.concatenate((Rp,Vp))
    return Xp

#
# Integra a equação de dinâmica perturbada por terceiro corpo para determinar a posição e
# velocidade heliocêntrica do VE a partir de agora é 100 dias solares médios
X0=np.concatenate((R0,V0))
T=100*24*60*60 # s - tempo final
sol=solve_ivp(dinPert3corpo,[0, T],X0,max_step=24*60*60)
#

# Propagação da órbita do VE usando a solução analítica kepleriana. Sol como primário.
# Ignorando a influência da Terra (terceiro corpo)
# Identificação dos parâmetros orbitais


par_orb=calculos_orbitais.determina_parametros_orbitais(0, mu, R0, V0,)
N=200
t=np.linspace(0, T, N)
Rk=np.empty([3,N]);Vk=np.empty([3,N])
for i in range(N):
    _,Rk[:,i],Vk[:,i]=propagacao.propagaEliptica(t[i],mu,par_orb)







plt.close("all")
#
plt.figure(1);
#
plt.subplot(2,3,1);plt.plot(sol.t/(24*60*60),sol.y[0]/1E9,label="integração")
plt.plot(t/(24*60*60),Rk[0,:]/1E9,label="kepleriana")
plt.xlabel('t(dms)');plt.ylabel('R_X (10^6 km)');plt.grid();plt.legend()
#
plt.subplot(2,3,2);plt.plot(sol.t/(24*60*60),sol.y[1]/1E9,label="integração")
plt.plot(t/(24*60*60),Rk[1,:]/1E9,label="kepleriana")
plt.xlabel('t(dms)');plt.ylabel('R_Y (10^6 km)');plt.grid();plt.legend()
#
plt.subplot(2,3,3);plt.plot(sol.t/(24*60*60),sol.y[2]/1E9,label="integração")
plt.plot(t/(24*60*60),Rk[2,:]/1E9,label="kepleriana");plt.legend()
plt.xlabel('t(dms)');plt.ylabel('R_Z (10^6 km)');plt.grid()
#
plt.subplot(2,3,4);plt.plot(sol.t/(24*60*60),sol.y[3]/1E3,label="integração")
plt.plot(t/(24*60*60),Vk[0,:]/1E3,label="kepleriana");plt.legend()
plt.xlabel('t(dms)');plt.ylabel('V_X (km/s)');plt.grid()
#
plt.subplot(2,3,5);plt.plot(sol.t/(24*60*60),sol.y[4]/1E3,label="integração")
plt.plot(t/(24*60*60),Vk[1,:]/1E3,label="kepleriana")
plt.xlabel('t(dms)');plt.ylabel('V_Y (km/s)');plt.grid();plt.legend()
#
plt.subplot(2,3,6);plt.plot(sol.t/(24*60*60),sol.y[5]/1E3,label="integração")
plt.plot(t/(24*60*60),Vk[2,:]/1E3,label="kepleriana")
plt.xlabel('t(dms)');plt.ylabel('V_Z (km/s)');plt.grid();plt.legend()
plt.show()