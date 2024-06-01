import numpy as np

def Vrel2Vine(vr, phir, Ar, we, r, dt):
    """
    Função para converter a velocidade, elevação e azimute da velocidade
    relativa para a velocidade inercial

    Entradas:
    vr (m/s): velocidade relativa
    phir (rad): inclinação da velocidade relativa
    Ar (rad): azimute da velocidade relativa
    we (rad/s): velocidade de rotação do referencial girante
    r (m): distância radial até a origem do referencial inercial
    dt (rad): latitude

    Saídas:
    vi (m/s): magnitude da velocidade com respeito ao referencial inercial
    phii (rad): ângulo de elevação da velocidade inercial (angulo de trajetoria)
    Ai (rad): ângulo de azimute da velocidade inercial
    """

    # Cálculos
    Ai = np.arctan2(vr * np.cos(phir) * np.sin(Ar) + we * r * np.cos(dt), vr * np.cos(phir) * np.cos(Ar))
    Ai = Ai if Ai >= 0 else Ai + 2 * np.pi

    vi = np.sqrt(vr*2 + 2 * vr * np.cos(phir) * np.sin(Ar) * r * we * np.cos(dt) + r**2 * we*2 * np.cos(dt)**2)

    phii = np.arctan2(np.sin(phir) * np.cos(Ai), np.cos(phir) * np.cos(Ar))
    if abs(phii) > np.pi/2:
        phii = phii + np.pi if phii < np.pi/2 else phii - np.pi

    return vi, phii, Ai