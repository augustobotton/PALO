import numpy as np
import math

def det_orbita(t0, rc0, vc0, mu):
    """
    Função para determinar parâmetros orbitais a partir de uma observação de
    posição e outra de velocidade, sendo as mesmas tomadas em relação ao
    primário em um problema de dois corpos e escritas no referencial
    celestial.
    """
    # Cálculos
    # Distância radial ao primário no instante observado
    r0 = np.linalg.norm(rc0)
    
    # Vetor quantidade de movimento angular específica no referencial celeste
    hc = np.cross(rc0, vc0)
    
    # Vetor excentricidade no sistema celeste
    ec = np.cross(vc0, hc) / mu - rc0 / r0
    
    # Excentricidade da órbita
    e = np.linalg.norm(ec)
    
    # Módulo do vetor hc
    h = np.linalg.norm(hc)
    
    # Parâmetro da órbita dada
    p = h**2 / mu
    
    # Semi eixo maior
    a = p / (1 - e**2)
    
    # Vetor parâmetro no referencial celeste
    pc = p * np.cross(hc, ec) / (h * e)
    
    # Anomalia verdadeira
    costheta = (p - r0) / (e * r0)
    sintheta = np.dot(rc0, pc) / (r0 * p)
    theta = math.atan2(sintheta, costheta)
    
    # O tempo de perigeu depende do tipo de órbita
    if (0 <= e) and (e < 1):
        tipo = 'e'  # órbita elíptica
    elif e == 1:
        tipo = 'p'  # órbita parabólica
    else:
        tipo = 'h'  # órbita hiperbólica
    
    # Tempo de perigeu
    if tipo == 'e':  # órbita elíptica
        # Movimento médio
        n = math.sqrt(mu / a**3)
        # Anomalia excêntrica
        E = 2 * math.atan(math.sqrt((1 - e) / (1 + e)) * math.tan(theta / 2))
        tau = t0 - (E - e * math.sin(E)) / n
    elif tipo == 'p':  # órbita parabólica
        tau = -((math.tan(theta / 2))*3 + 3 * math.tan(theta / 2)) / (mu / p*3)**(1 / 6)
    else:  # órbita hiperbólica
        # Movimento médio hiperbólico
        n = math.sqrt(-mu / a**3)
        # Anomalia hiperbólica
        H = 2 * math.atanh(math.sqrt((e - 1) / (1 + e)) * math.tan(theta / 2))
        tau = -(e * math.sinh(H) - H) / n
    
    # Linha dos nodos
    # Vetor unitário ao longo do vetor h (no sistema celeste)
    ih = hc / h
    # Vetor unitário ao longo da linha dos nodos (no sistema celeste)
    Kc = np.array([0, 0, 1])
    nc = np.cross(Kc, ih) / np.linalg.norm(np.cross(Kc, ih))
    # Ascensão reta do nodo ascendente
    OMEGA = math.atan2(nc[2], nc[1])
    # Inclinação
    i = math.acos(np.dot(ih, Kc))
    # Vetor unitário ao longo do vetor excentricidade (no referencial celeste)
    ie = ec / e
    # Argumento de perigeu
    cosomega = np.dot(ie, nc)
    sinomega = np.dot(ih, np.cross(nc, ie))
    omega = math.atan2(sinomega, cosomega)
    
    # Vetor de parâmetros de saída
    par_orb = np.array([a, e, tau, OMEGA, i, omega])
    return par_orb