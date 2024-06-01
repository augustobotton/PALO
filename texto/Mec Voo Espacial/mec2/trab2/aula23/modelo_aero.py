import numpy as np
from scipy import interpolate
# Modelo de arrasto em função dos números de Mach e Knudsen sem levar em conta
# os ângulos de ataque e derrapagem
#
# Função para interpolar o coeficiente de arrasto do meio contínuo em função do
# número de Mach
def iterpCD(M):
# Função para calcular o coeficiente de arrasto contínuo
# Dados extraídos de Tewari
    CD=np.array([[4.44e-16, 0.47],
    [0.16, 0.47],[0.34, 0.48],[0.48, 0.49],[0.53, 0.50],[0.64, 0.51],[0.71, 0.54],
    [0.75, 0.56],[0.85, 0.59],[0.87, 0.62],[0.91, 0.67],[0.94, 0.64],[0.95, 0.70],
    [0.95, 0.72],[0.97, 0.79],[0.98, 0.80],[0.99, 0.87],[1.0, 0.89],[1.01, 1.03],
    [1.011, 0.98],[1.018, 0.95],[1.025, 1.013],[1.034, 0.97],[1.13, 1.013],
    [1.15, 1.0],[1.22, 0.98],[1.26, 0.97],[1.27, 0.96],[1.30, 0.94],[1.35, 0.92],
    [1.40, 0.90],[1.44, 0.88],[1.49, 0.86],[1.53, 0.84],[1.58, 0.82],[1.67, 0.79],
    [1.74, 0.77],[1.81, 0.75],[1.88, 0.73],[1.98, 0.70],[2.18, 0.66],[2.28, 0.64],
    [2.36, 0.63],[2.61, 0.59],[2.73, 0.58],[2.87, 0.56],[2.99, 0.55],[3.17, 0.54],
    [3.43, 0.53],[3.68, 0.52],[3.97, 0.51],[4.18, 0.51],[4.44, 0.50],[4.69, 0.50],
    [4.95, 0.49],[5.18, 0.49],[5.44, 0.49],[5.69, 0.49],[5.95, 0.49],[6.20, 0.49],
    [6.45, 0.49],[6.70, 0.49],[6.95, 0.49],[7.20, 0.49],[7.46, 0.49],[7.71, 0.49],
    [7.96, 0.49],[8.21, 0.49],[8.46, 0.49],[8.71, 0.49],[8.97, 0.48],[9.22, 0.48],
    [9.47, 0.48],[9.72, 0.48],[10.22, 0.48],[10.47, 0.48],[10.73, 0.48],
    [10.98, 0.48],[11.23, 0.48],[11.48, 0.48],[11.73, 0.48]]);
    # Interpolação para o número de Mach da entrada
    f=interpolate.interp1d(CD[:,0],CD[:,1])
    CDc=f(M)
    return CDc

#
def coefArrastoFoguete(V,h,M,Kn,T,R):
    # Modelo de arrasto conforme a referencia
    # TEWARI, A. Atmospheric and Space Flight Dynamics:
    # Modelling and simulation with MATLAB and Simulink. Boston: Birkhauser, 2007.
    # Exemplo 12.6
    # Entradas
    # V: velocidade (m/s)
    # M: Número de Mach
    # Kn: Número de Knudsen
    # T: Temperatura (K)
    # R: Constante de gas ideal do ar (J/kg.K)
    # Saída
    # CD: Coeficiente de arrasto
    s=V/np.sqrt(2*R*T);
    CDfm=1.75+np.sqrt(np.pi)/(2*s);
    CDc=iterpCD(M);
    if h<2000e3:
        if Kn<0.0146: # Meio contínuo
            CD=CDc;
        elif Kn<14.5: # Transição entre o contínuo e o livre molecular
            CD=CDc+(CDfm-CDc)*((1/3)*np.log10(Kn/np.sin(30*np.pi/180))+0.05113);
        else: # Escoamento livre molecular
            CD=CDfm
    else: # Sem arrasto
        CD=0;
    return CD