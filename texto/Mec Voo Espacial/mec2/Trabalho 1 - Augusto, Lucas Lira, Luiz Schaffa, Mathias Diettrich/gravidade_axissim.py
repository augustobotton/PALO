#Adaptado de André Luis da Silva
import numpy as np 
import parametros 


def grav_axissimetrico(r,phi):
# Entradas 
  
#r: Distância Radial ao centro do massa do Planeta(m).
# phi: Colatitude (rad)
# 
# Saídas:
# gr: componente radial da gravidade [m/s^2].
# gphi: componente colatitudinal (sul) da gravidade[m/s^2] 

# Paramêtros do modelo:
# 
# 
    Re = parametros.RE
    G =  parametros.G
    M  = parametros.M
    J2 = parametros.J2; J3=parametros.J3; J4=parametros.J4

    gr =  G*M*(2*J2*Re**2*(1.5*np.cos(phi)**2 - 0.5)/r**3 + 3*J3*Re**3*(2.5*np.cos(phi)**3 - 1.5*np.cos(phi))/r**4 + 4*J4*Re**4*(4.375*np.cos(phi)**4 - 3.75*np.cos(phi)**2 + 0.375)/r**5)/r - G*M*(-J2*Re**2*(1.5*np.cos(phi)**2 - 0.5)/r**2 - J3*Re**3*(2.5*np.cos(phi)**3 - 1.5*np.cos(phi))/r**3 - J4*Re**4*(4.375*np.cos(phi)**4 - 3.75*np.cos(phi)**2 + 0.375)/r**4 + 1)/r**2
    
    gphi =  G*M*(3.0*J2*Re**2*np.sin(phi)*np.cos(phi)/r**2 - J3*Re**3*(-7.5*np.sin(phi)*np.cos(phi)**2 + 1.5*np.sin(phi))/r**3 - J4*Re**4*(-17.5*np.sin(phi)*np.cos(phi)**3 + 7.5*np.sin(phi)*np.cos(phi))/r**4)/r**2

    return gr, gphi