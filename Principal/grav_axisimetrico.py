import numpy as np
import parametros

def grav_axisimetrico(r, delta):
    # Entradas:
    # r: distância radial ao centro de massa do planeta (m).
    # delta: latitude centrada no planeta (rad)
    # Saidas:
    # gc: componente centripeta da gravidade (m/s^2).
    # gdel: componente latitudinal (norte) da gravidade (m/s^2)
    ## Variaveis
    Re = parametros.Re
    mut = parametros.mut
    J2 = parametros.J2
    J3 = parametros.J3
    J4 = parametros.J4

    # Colatitude (rad)
    phi = np.pi/2 - delta
    # Componente radial da gravidade de corpo axisimetrico
    gr = -(mut) / r**2 - (1/r) * mut * (-((J2*Re**2*(-1 + 3*np.cos(phi)**2))/r**3) - (1.5*J3*Re**3*(-3*np.cos(phi) + 5*np.cos(phi)**3))/r**4 -
          (J4*Re**4*(3-30*np.cos(phi)**2 + 35*np.cos(phi)**4))/(2*r**5)) + (mut/r**2) * ((0.5*J2*Re**2*(-1+3*np.cos(phi)**2))/r**2 + (0.5*J3*Re**3*(-3*np.cos(phi) + 5*np.cos(phi)**3))/r**3 +
          (J4*Re**4*(3-30*np.cos(phi)**2+35*np.cos(phi)**4))/(8*r**4))
    # Componente sul da gravidade de corpo axisimetrico
    gphi = -(1/(r**2)) * mut * (-((3*J2*Re**2*np.cos(phi)*np.sin(phi))/r**2) + (0.5*J3*Re**3*(3*np.sin(phi)-15*np.cos(phi)**2*np.sin(phi)))/r**3 +
            (J4*Re**4*(60*np.cos(phi)*np.sin(phi)-140*np.cos(phi)**3*np.sin(phi)))/(8*r**4))
    # Componentes nas direcoes centripeta e latitudinal (norte)
    gc = -gr
    gdel = -gphi
    return gc, gdel
