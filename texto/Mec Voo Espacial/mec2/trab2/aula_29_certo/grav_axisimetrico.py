from principal_GSO import *

def grav_axisimetrico(r,delta):

    phi=np.pi/2-delta
    
    gr = -(mut) / (r ** 2) - (1 / r) * mut * (-((J2 * Requat ** 2 * (-1 + 3 * np.cos(phi) ** 2)) / (r ** 3)) - (1.5 * J3 * Requat ** 3 * (-3 * np.cos(phi) +
    5 * np.cos(phi) ** 3)) / (r ** 4) -(J4 * Requat ** 4 * (3 - 30 * np.cos(phi) ** 2 + 35 * np.cos(phi) ** 4)) / (2 * r ** 5)) + \
    (mut / (r ** 2)) * ((0.5 * J2 * Requat ** 2 * (-1 + 3 * np.cos(phi) ** 2)) / (r ** 2) + (0.5 * J3 * Requat ** 3 * (-3 * np.cos(phi) + \
    5 * np.cos(phi) ** 3)) / (r ** 3) + (J4 * Requat ** 4 * (3 - 30 * np.cos(phi) ** 2 + 35 * np.cos(phi) ** 4)) / (8 * r ** 4))

    gphi = -(1. / np.square(r)) * mut * (
    -(3. * J2 * np.square(Requat) * np.cos(phi) * np.sin(phi)) / np.power(r, 2)
    + (0.5 * J3 * np.power(Requat, 3) * (3. * np.sin(phi) - 15. * np.square(np.cos(phi)) * np.sin(phi))) / np.power(r, 3)
    + (J4 * np.power(Requat, 4) * (60. * np.cos(phi) * np.sin(phi) - 140. * np.power(np.cos(phi), 3) * np.sin(phi))) / (8. * np.power(r, 4)))

    gc = -gr
    gdel=-gphi

    return gc, gdel
 