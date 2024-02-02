import numpy as np


class Converte:

    def long_ECEF2ECI(self, t, long, we, tg):
        # Função para calcular a longitude celeste a partir da longitude fixa ao planeta
        # Entradas
        # t (s) - Tempo no qual se deseja saber a longitude celeste
        # long (rad) - Longitude relativa ao referencial fixo ao planeta
        # we (rad/s) - Velocidade de rotação do planeta
        # tg (s) - Tempo no qual o meridiano de referência tem longitude celeste nula
        # Saída
        # long_c (rad) - Longitude celeste no tempo t
        long_c = long + we * (t - tg)
        return long_c

    def RvelPolar2RvelRet(self, v, A, phi, r, lat, lon):
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

        V = CLH.T @ Vlvlh[:][:]

        # Vetor posição no sistema LHVLH
        Rlvlh = np.array([[r],
                          [0],
                          [0]])
        # Transformação da posição para o sistema ECI ou ECEF em coordenadas
        # retangulares
        R = CLH.T @ Rlvlh
        return R, V

    def Vrel2Vine(self, vr, phir, Ar, we, r, dt):
        # Função para converter a velocidade, elevação e azimute da velocidade relativa para a velocidade inercial
        # Entradas
        # vr (m/s): velocidade relativa
        # phir (rad): inclinação da velocidade relativa
        # Ar (rad): azimute da velocidade relativa
        # we (rad/s): velocidade de rotação do referencial girante
        # r (m): distância radial até a origem do referencial inercial
        # dt (rad): latitude
        # Saídas
        # v (m/s): magnitude da velocidade com respeito ao referencial inercial
        # phi (rad): ângulo de elevação da velocidade inercial (angulo de trajetoria)
        # A (rad): ângulo de azimute da velocidade inercial
        ## Cálculos
        Ai = np.arctan2(vr * np.cos(phir) * np.sin(Ar) + we * r * np.cos(dt), vr * np.cos(phir) * np.cos(Ar))
        if Ai < 0:
            Ai += 2 * np.pi
        vi = np.sqrt(
            vr ** 2 + 2 * vr * np.cos(phir) * np.sin(Ar) * r * we * np.cos(dt) + r ** 2 * we ** 2 * np.cos(dt) ** 2)
        phii = np.arctan((np.sin(phir) * np.cos(Ai)) / (np.cos(phir) * np.cos(Ar)))
        return vi, phii, Ai
