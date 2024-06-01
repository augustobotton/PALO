import numpy as np

def RvelPolar2RvelRet(v, A, phi, r, lat, long):
    """
    Função para converter velocidade do sistema LVLH (coordenadas polares)
    para o sistema ECI ou ECEF retangular
    A velocidade pode ser a relativa ou a inercial, o resultado final será
    correlato

    Entradas:
    v (m/s): Módulo do vetor velocidade
    A (rad): Azimute da velocidade
    phi (rad): Elevação da velocidade
    r (m): distância radial
    lat (rad): Latitude 
    long (rad): Longitude no referencial desejado (ECI ou ECEF)

    Saídas:
    R (m): Vetor posição em coordenadas retangulares no
    sistema ECI ou ECEF (dependendo da entrada de dados de longitude)
    V (m/s): Vetor velocidade em coordenadas retangulares no
    sistema ECI ou ECEF (dependendo da entrada de dados de longitude). Pode
    ser a velocidade relativa ou a inercial, dependendo dos dados de
    velocidade fornecidos.
    """
    # Matriz de conversão do sistema ECI ou ECEF para o LVLH 
    CLH = np.array([[np.cos(lat)*np.cos(long), np.cos(lat)*np.sin(long), np.sin(lat)],
                    [-np.sin(long), np.cos(long), 0],
                    [-np.sin(lat)*np.cos(long), -np.sin(lat)*np.sin(long), np.cos(lat)]])

    # Vetor velocidade em coordenadas cartezianas no sistema LVLH
    Vlvlh = v * np.array([np.sin(phi),
                          np.cos(phi)*np.sin(A),
                          np.cos(phi)*np.cos(A)])

    # Transformação da velocidade para o sistema ECI ou ECEF em coordenadas retangulares
    V = np.dot(CLH.T, Vlvlh)

    # Vetor posição no sistema LHVLH
    Rlvlh = np.array([r, 0, 0])

    # Transformação da posição para o sistema ECI ou ECEF em coordenadas retangulares
    R = np.dot(CLH.T, Rlvlh)

    return R,V