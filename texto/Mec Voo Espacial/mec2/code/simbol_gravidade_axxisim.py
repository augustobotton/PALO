#Script para calcular a gravidade de um planeta Axissimétrico com 3 cnst de Jeffrey
#Cálculo Simbólico - Adaptado de André Luis da Silva

import sympy as sym
#Definição das funções para cálculo dos polinômios de Lagrange
def P2(x):
    return (1/2)*(3*x**2-1)

def P3(x): 
    return (1/2)*(5*x**3-3*x)

def P4(x): 
    return (1/8)*(35*x**4-30*x**2+3)


#Declara Símbolos
r, phi, G, M, Re, J2, J3, J4 = sym.symbols('r phi G M Re J2 J3 J4')

#Potencial Gravitacional 
PHI = (G*M/r)*(1-((Re/r)**2*J2*P2(sym.cos(phi))+(Re/r)**3*J3*P3(sym.cos(phi))+(Re/r)**4*J4*P4(sym.cos(phi))))
print('PHI =', PHI)

#Cálculo da Aceleração da Gravidade em Coordenadas Esféricas
gr = sym.diff(PHI,r)
gphi = sym.diff(PHI,phi)/r
print('gr = ', gr)
print('gphi = ', gphi)
