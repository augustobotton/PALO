import numpy as np
from scipy.optimize import fsolve
import sistReferencia
 #
 # Módulo com funções de propagação de órbita kepleriana
 #
 # Matriz de transição de estado
def matrizTransicaoEstado(theta,theta0,mu,orb):
 # Funcao para determinar a matriz de transicao de estado cujos elementos
 # sao os coeficientes de Lagrange
 # Entradas:
 # theta [rad]: anomalia verdadeira
 # mu: [m^3/s^2] parâmetro gravitacional do primário.
 # orb: vetor 6x1 de elementos orbitais clássicos
 # orb[0]: a [m] - semi eixo maior
 # orb[1]: e - excentricidade
 # orb[2]: tau [s] - tempo de periastro
 # Saída:
 # PHI [2x2]: matriz de transicao de estado
 # Elementos orbitais
 a=orb[0];e=orb[1];tau=orb[2]
 # Calculo das constantes associadas a orbita
 p=(1-e**2)*a # [m] Parametro
 h=np.sqrt(p*mu) # momento angular especifico da orbita
 r0=p/(1+e*np.cos(theta0)) # distancia radial inicial
 # Calculo dos coeficientes de Lagrange
 r=p/(1+e*np.cos(theta)) # Distancia radial para o theta dado
 f=1+(r/p)*(np.cos(theta-theta0)-1)
 g=(r*r0/h)*np.sin(theta-theta0)
 fp=-(h/p**2)*(np.sin(theta-theta0)+e*(np.sin(theta)-np.sin(theta0)))
 gp=1+(r0/p)*(np.cos(theta-theta0)-1)
 # Matriz de transicao de estado
 PHI=np.array([[f,g],[fp,gp]])
 return PHI
 #

 # Órbita elíptica: cálculo da propagação temporal da anomalia excêntrica
 #
def Kepler(E,*dados):
 # Funcao objetivo: Equação de Kepler
 # Entradas:
 # E: [rad] anomalia excêntrica - incógnita
 # *dados: parâmetros da função
 # Saída:
 # y: quando igual a zero, a equacao de Kepler esta resolvida
 #
 # Recebimento de parametros
 e,M=dados
 # Quando y=0, a equacao de Kepler esta resolvida
 y=E-e*np.sin(E)-M
 return y
 
 #

def resolveEqKepler(t,mu,orb):
 # Órbita eliptica: calcula a anomalia excentrica para dado tempo
 # Entrada
 # t: [s] tempo
 # mu: [m^3/s^2] parâmetro gravitacional do primário
 # orb: vetor de elementos orbitais clássicos
 # orb(0): a [m] - semi eixo maior
 # orb(1): e - excentricidade
 # orb(2): tau [s] - tempo de periastro, medido com respeito a t=0
 # Saida
 # E: [rad] anomalia excentrica no instante t
 # Elementos orbitais
 a = orb[0] # Semi eixo maior
 e = orb[1] # Excentricidade
 tau = orb[2] # tempo de periastro
 # Movimento medio
 n = np.sqrt(mu/a**3)
 # Anomalia media
 M = n*(t-tau)
 # Passagem de parâmetros para a função objetivo
 dados = (e,M)
 # Chute inicial
 E0 = M # Seria o resultado da orbita circular
 # Resolve numericamente a equacao de Kepler com a funcao fsolve do pacote
 # scipy.optimize
 E = fsolve(Kepler,E0,args=dados)
 return E[0]
 # fsolve é bandida, temos que indexar para obter a solução em float, não ela
 # retorna uma tupla, mesmo que a solução seja só um valor.
 #
 # Função para propagação de órbita elíptica no referencial inercial
def propagaEliptica(t,mu,orb):
 # Entradas:
 # t: [s] tempo. A sua referência é o tempo de periastro fornecido.
 # mu: [m^3/s^2] parâmetro gravitacional do primário.
 # orb: vetor 6x1 de elementos orbitais clássicos
 # orb[0]: a [m] - semi eixo maior
 # orb[1]: e - excentricidade
 # orb[2]: tau [s] - tempo de periastro, medido em relação ao tempo t=0
 # orb[3]: OMEGA [rad] - longitude celeste do nodo ascendente do referencial
 # perifocal com respeito ao inercial.
 # orb[4]: i [rad] - inclinação da órbita com respeito ao plano XY do referencialinercial
 # orb[5]: omega [rad] - argumento de periastro
 # Saídas:
 # Variáveis calculadas no instante t fornecido
 # theta: [rad] - anomalia verdadeira
 # Ri: [m] - Vetor 3x1, posição no referencial inercial. Coordenadas retangulares
 # Vi: [m] - Vetor 3x1, velocidade no referencial inercial. Coordenadas retangulares
 # Elementos orbitais
 a=orb[0];e=orb[1];tau=orb[2];OMEGA=orb[3];inc=orb[4];omega=orb[5]
 # Anomalia verdadeira em t=0
 E0=resolveEqKepler(0,mu,orb)
 theta0=2*np.arctan(np.sqrt((1+e)/(1-e))*np.tan(E0/2))
 # Parametro
 p=(1-e**2)*a
 # Período
 n=np.sqrt(mu/a**3);f=n/(2*np.pi);P=1/f
 # Vetores posição e velocidade inicial, no referencial perifocal, escritos
 # em coordenadas retangulares
 h=np.sqrt(p*mu) # [m^2/s] Quantidade de movimento angular especifica
 r0=p/(1+e*np.cos(theta0))
 R0=np.array([r0*np.cos(theta0),r0*np.sin(theta0),0])
 V0=np.array([-(mu/h)*np.sin(theta0),(mu/h)*(e+np.cos(theta0)),0])
 # Resolve a equacao de Kepler, determinando a anomalia excêntrica
 E=resolveEqKepler(t,mu,orb)
 # Determina a anomalia verdadeira
 theta=2*np.arctan(np.sqrt((1+e)/(1-e))*np.tan(E/2))
 if theta<0:
    theta=theta+2*np.pi
 # Determina a matriz de transição de estado a partir de theta
 PHI=matrizTransicaoEstado(theta,theta0,mu,orb)
 # Determina posição e velocidade perifocal em função da anomalia
 # verdadeira pela matriz de transição de estado
 R=PHI[0,0]*R0+PHI[0,1]*V0
 V=PHI[1,0]*R0+PHI[1,1]*V0
 # Matriz de transformação de coordenadas do referencial inercial para o perifocal
 Cip=sistReferencia.matInercPerif(OMEGA,inc,omega)
 # Do perifocal para o inercial
 Cpi=np.transpose(Cip)
 # Posição e velocidade no referencial inercial
 Ri=Cpi@R;Vi=Cpi@V
 return theta,Ri,Vi

mu=1.32712440018e20
at=149597870E3 # m
et= 0.01667
taut=-100*24*60*60 # s
orbt=np.array([at,et,taut,0,0,0])
propagaEliptica(0,mu,orbt)