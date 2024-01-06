import numpy as np

def RvelPolar2RvelRet(v, A, phi, r, lat, lon):
    """
    Função para converter velocidade do sistema LVLH (coordenadas polares) para o sistema ECI ou ECEF retangular
    """
    # Entradas:
    # v (m/s): Módulo do vetor velocidade
    # A (rad): Azimute da velocidade
    # phi (rad): Elevação da velocidade
    # r (m): distância radial
    # lat (rad): Latitude
    # lon (rad): Longitude no referencial desejado (ECI ou ECEF)
    # Saídas:
    # R=[R_X, R_Y, R_Z]^T (m): Vetor posição em coordenadas retangulares no sistema ECI ou ECEF (dependendo da entrada de dados de longitude)
    # V=[V_X, V_Y, V_Z]^T (m/s): Vetor velocidade em coordenadas retangulares no sistema ECI ou ECEF (dependendo da entrada de dados de longitude). Pode ser a velocidade relativa ou a inercial, dependendo dos dados de velocidade fornecidos.

    ## Cálculos
    # Matriz de conversão do sistema ECI ou ECEF para o LVLH
    CLH = np.array([[np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)],
                    [-np.sin(lon), np.cos(lon), 0],
                    [-np.sin(lat) * np.cos(lon), -np.sin(lat) * np.sin(lon), np.cos(lat)]], dtype=('object'));
    # Vetor velocidade em coordenadas cartezianas no sistema LVLH
    Vlvlh = v * np.array([[np.sin(phi)],
                          [np.cos(phi) * np.sin(A)],
                          [np.cos(phi) * np.cos(A)]], dtype=('object'));
    # Transformação da velocidade para o sistema ECI ou ECEF em coordenadas retangulares

    V = CLH.T@Vlvlh[:][:]

    # Vetor posição no sistema LHVLH
    Rlvlh = np.array([[r],
                      [0],
                      [0]])
    # Transformação da posição para o sistema ECI ou ECEF em coordenadas
    # retangulares
    R = CLH.T@Rlvlh
    return R, V
