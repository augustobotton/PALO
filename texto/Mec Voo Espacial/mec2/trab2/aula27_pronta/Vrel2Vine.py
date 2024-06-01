import numpy as np

def Vrel2Vine(vr, phir, Ar, we, r, dt):
    # Função para converter a velocidade, elevação e azimute da velocidade
    # relativa para a velocidade inercial
    
    # Cálculos
    Ai = np.arctan2(vr*np.cos(phir)*np.sin(Ar) + we*r*np.cos(dt), vr*np.cos(phir)*np.cos(Ar))
    if Ai < 0:
        Ai = Ai + 2*np.pi

    vi = np.sqrt((vr**2) + 2*vr*np.cos(phir)*np.sin(Ar)*r*we*np.cos(dt) + (r**2)*(we**2)*(np.cos(dt)**2))
    phii = np.arctan2(np.sin(phir)*np.cos(Ai), np.cos(phir)*np.cos(Ar))
    
    return vi, phii, Ai