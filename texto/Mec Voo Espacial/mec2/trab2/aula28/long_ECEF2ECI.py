def long_ECEF2ECI(t, long, we, tg):
    # Função para calcular a longitude celeste a partir da longitude
    # fixa ao planeta
    
    # Cálculo
    long_c = long + we * (t - tg)
    return long_c