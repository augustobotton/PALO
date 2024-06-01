import numpy as np

def RvelPolar2RvelRet(v, A, phi, r, lat, long):
    # Funcao para converter velocidade do sistema LVLH (coordenadas polares)
    # para o sistema ECI ou ECEF retangular
    # A velocidade pode ser a relativa ou a inercial, o resultado final sera
    # correlato

    # Calculos
    # Matriz de conversao do sistema ECI ou ECEF para o LVLH 
    CLH = np.array([[np.cos(lat)*np.cos(long), np.cos(lat)*np.sin(long), np.sin(lat)],
                    [-np.sin(long), np.cos(long), 0],
                    [-np.sin(lat)*np.cos(long), -np.sin(lat)*np.sin(long), np.cos(lat)]])

    # Vetor velocidade em coordenadas cartesianas no sistema LVLH
    Vlvlh = v * np.array([np.sin(phi), np.cos(phi)*np.sin(A), np.cos(phi)*np.cos(A)])

    # Transformacao da velocidade para o sistema ECI ou ECEF em coordenadas retangulares
    V = np.dot(CLH.T, Vlvlh)

    # Vetor posicao no sistema LHVLH
    Rlvlh = np.array([r, 0, 0])

    # Transformacao da posicao para o sistema ECI ou ECEF em coordenadas
    # retangulares
    R = np.dot(CLH.T, Rlvlh)

    return R, V

