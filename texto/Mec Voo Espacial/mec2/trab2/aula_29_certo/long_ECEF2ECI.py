import numpy as np

def long_ECEF2ECI(t, long, we, tg):
    """
    Função para calcular a longitude celeste a partir da longitude
    fixa ao planeta
    
    Parâmetros:
    t (s) - Tempo no qual se deseja saber a longitude celeste
    long (rad) - Longitude relativa ao referencial fixo ao planeta
    we (rad/s) - Velocidade de rotação do planeta
    tg (s) - Tempo no qual o meridiano de referência tem longitude celeste
    nula
    
    Retorno:
    long_c (rad) - Longitude celeste no tempo t
    """
    long_c = long + we * (t - tg)
    return long_c