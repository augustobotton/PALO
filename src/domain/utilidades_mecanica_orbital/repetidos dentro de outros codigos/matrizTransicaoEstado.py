import parametros
import numpy as np
def matrizTransicaoEstado(theta):
# Funcao para determinar a matriz de transicao de estado cujos elementos
# sao os coeficientes de Lagrange
# Entrada
# theta [rad]: anomalia verdadeira
# Saida
# PHI [2x2]: matriz de transicao de estado
#
# As variaveis globais sao passadas por um arquivo de parametros
    theta0=parametros.THETA0
    e=parametros.E
    p=parametros.P
    mu=parametros.MU
# Calculo das constantes associadas a orbita
    h=np.sqrt(p*mu) # momento angular especifico da orbita
    r0=p/(1+e*np.cos(theta0))
# distancia radial inicial
# Calculo dos coeficientes de Lagrange
    r=p/(1+e*np.cos(theta))
# Distancia radial para o theta dado
    f=1+(r/p)*(np.cos(theta-theta0)-1)
    g=(r*r0/h)*np.sin(theta-theta0)
    fp=-(h/p**2)*(np.sin(theta-theta0)+e*(np.sin(theta)-np.sin(theta0)))
    gp=1+(r0/p)*(np.cos(theta-theta0)-1)
# Matriz de transicao de estado
    PHI=np.array([[f,g],[fp,gp]])
    return PHI
