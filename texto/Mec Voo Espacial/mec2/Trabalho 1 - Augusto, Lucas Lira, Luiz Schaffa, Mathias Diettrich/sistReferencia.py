# Pacote com funções para cálculos associados a sistemas de referência
import numpy as np
#
# Matrizes de rotação elementares
#
def rotX(alf):
    C1=np.array([[1,0,0],
    [0,np.cos(alf),np.sin(alf)],
    [0,-np.sin(alf),np.cos(alf)]])
    return C1
#

def rotY(alf):
    C2=np.array([[np.cos(alf),0,-np.sin(alf)],
    [0,1,0],
    [np.sin(alf),0,np.cos(alf)]])
    return C2
#
def rotZ(alf):
    C3=np.array([[np.cos(alf),np.sin(alf),0],
    [-np.sin(alf),np.cos(alf),0],
    [0,0,1]])
    return C3
#
# Rotação do sistema de refrência inercial para o perifocal
def matInercPerif(OMEGA,inc,omega):
    Cip=rotZ(omega)@rotX(inc)@rotZ(OMEGA)
    return Cip