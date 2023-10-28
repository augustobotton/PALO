import numpy as np


def atm_padrao(h, v, lc, dT):
    """
    Funcao para calculo da atmosfera padrao para a altitude geometrica de zero ate 2000 km
    """
    # Entradas:
    # h - altitude geometrica [m]
    # v - velocidade do veiculo em relacao ao escoamento nao perturbado, em metros por segundo [m/s]
    # lc - comprimento caracteristico do veiculo, em metros [m]
    # dT - Variacao de temperatura em relacao a temperatura padrao [K ou °C]
    # Saidas:
    # T - Temperatura, em Kelvin [K]
    # p - Pressao, em Pascal [N/m^2]
    # rho - Densidade [kg/m^2]
    # ainf - Velocidade do som [m/s]
    # M - Numero de Mach [adm]
    # mu - Coeficiente de viscosidade dinamica [kg/m*s]
    # Pr - Numero de Prandtl [adm]
    # Kn - Numero de Knudsen [adm]
    # d - Parametro de regime de escoamento [adm]
    # Re - Numero de Reynolds [adm]

    ## Entrada de dados
    # Vetores com os seguintes dados, cada linha eh uma camada do modelo de
    # atmosfera padrao: altitude no inicio da camada (m), temperatura no inicio
    # da camada (K), constante de gas ideal do ar na camada (J/kg.K), taxa de
    # lapso termico (K/m)

    hi = 1e3*np.array([0, 11.0191, 20.0631, 32.1619, 47.3501, 51.4125, 71.8020, 86, 100, 110, 120, 150, 160, 170, 190, 230, 300,
          400, 500, 600, 700]).T
    Ti = np.array([288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65, 186.946, 210.02, 257.0, 349.49, 892.79, 1022.2,
          1103.4, 1205.4, 1322.3, 1432.1, 1487.4, 1506.1, 1506.1, 1507.6]).T
    Ri = np.array([287.0, 287.0, 287.0, 287.0, 287.0, 287.0, 287.02, 287.02, 287.84, 291.06, 308.79, 311.80, 313.69, 321.57,
         336.68, 366.84, 416.88, 463.36, 493.63, 514.08, 514.08]).T
    a = 1e-3*np.array([-6.5, 0.0, 1.0, 2.8, 0.0, -2.8, -2.0, 1.693, 5.0, 10.0, 20.0, 15.0, 10.0, 7.0, 5.0, 4.0, 3.3, 2.6, 1.7, 1.1, 0.0]).T
    Mi = np.array([[28.9644], [28.9644], [28.9644], [28.9644], [28.9644], [28.9644], [28.9644], [28.9644],[28.88], [28.56], [28.07], [26.92], [26.66], [26.4], [25.85], [24.7], [22.66], [19.94], [17.94],[16.84], [16.17]])


    # Constantes
    M0 = Mi[0,0]
    g0 = 9.80665  # Valor ao nivel do mar da aceleracao da gravidade (m/s^2)
    Na = 6.0220978e23  # Numero de Avogadro
    sigma = 3.65e-10  # Diametro de colisao para o ar (m)
    m0 = 28.964e-3  # Massa molar do ar ao nivel do mar (kg/Mol)
    P0 = 1.01325e5  # Pressao padrao ao nivel do mar (N/m^2)
    Re = 6378.14e3  # Raio medio da Terra (m)
    gamma = 1.405  # Razao de calores especificos ao nivel do mar

    ## Constantes calculadas
    # Numero beta associado a distancia radial media do nivel do mar
    beta = 2 / Re

    ## Identifica a camada a qual a altitude pertence
    # Contador
    i = 0

    if h < 0:
        # print('Cuidado! Altitude negativa.')
        # print('Os resultados apresentados dizem respeito a h=0.')
        i = 0
        h = 0
    elif h > 2000e3:
        # print('A altitude fornecida esta acima do limite superior de 2.000 km.')
        # print('Os resultados apresentados dizem respeito a h=2.000 km')
        i = 20
        h = 2000e3
    else:
        for i in range(1, 22):
            if i == 21:
                break
            elif h >= hi[i - 1] and h < hi[i]:
                break

    ## Realiza os calculos
    # Pressao no inicio da camada i
    Pi = P0  # Pressao inicial da camada - inicializa com o valor ao nivel do mar
    px = P0  # Variavel auxiliar para guardar o valor de pressao inicial na camada anterior

    for j in range(2, i+1):
        if a[j - 2] != 0:  # Verifica se a camada nao eh isotermica
            # Calcula a pressao inicial da camada i a partir do modelo da camada i-1
            A = 1 + (a[j - 2] * (hi[j-1] - hi[j - 2])) / Ti[j - 2]
            B = -(g0 / (Ri[j - 1] * a[j - 2])) * (1 + beta * (Ti[j - 2] / a[j - 2] - hi[j - 2]))
            C = (g0 * beta / (Ri[j - 1] * a[j - 2])) * (hi[j-1] - hi[j - 2])
            Pi = px * (A ** B) * np.exp(C)
            px = Pi  # O valor atual sera o valor anterior na proxima iteracao
        else:
            # Calcula a pressao inicial da camada i pelo modelo isotermico
            Pi = px * np.exp(-(g0 / (Ri[j-1] * Ti[j - 2])) * (hi[j-1] - hi[j - 2]) * (1 - (beta / 2) * (hi[j-1] - hi[j - 2])))
            px = Pi  # O valor atual sera o valor anterior na proxima iteracao

    # Temperatura padrao
    Tm = Ti[i-1] + a[i-1] * (h - hi[i-1])

    # Corrige pelo valor de delta T adotado
    Tm = Tm + dT

    if i >= 21:
        R = Ri[i-1]
        M = Mi[i-1]
    else:
        R = Ri[i-1] + ((Ri[i] - Ri[i-1]) / (hi[i] - hi[i-1])) * (h - hi[i-1])
        M = Mi[i-1] + ((Mi[i] - Mi[i-1]) / (hi[i] - hi[i-1])) * (h - hi[i-1])

    # Temperatura padrao (cinetica) - saida da funcao
    T = (M / M0) * Tm

    # Pressao
    if a[i-1] != 0:  # Verifica se a camada não é isotérmica
        # Calcula a pressão
        A = 1 + (a[i-1] * (h - hi[i-1])) / Ti[i-1]
        B = -(g0 / (R * a[i-1])) * (1 + beta * (Ti[i-1] / a[i-1] - hi[i-1]))
        C = (g0 * beta / (R * a[i-1])) * (h - hi[i-1])
        p = Pi * (A ** B) * np.exp(C)
    else:
        # Calcula a pressão pelo modelo isotérmico
        p = Pi * np.exp(-(g0 / (R * Ti[i-1])) * (h - hi[i-1]) * (1 - (beta / 2) * (h - hi[i-1])))

    # Calcula a densidade
    rho = p / (R * Tm)

    # Velocidade do som
    ainf = np.sqrt(gamma * R * Tm)

    # Numero de Mach
    M = v / ainf

    # Coeficiente de viscosidade dinamica
    mu = 1.458e-6 * (Tm ** (3 / 2)) / (Tm + 110.4)

    # Numero de Prandtl
    cp = R * gamma / (gamma - 1)
    kT = (2.64638e-3 * (Tm ** (3 / 2))) / (Tm + 245.4 * (10 ** (-12 / Tm)))
    Pr = mu * cp / kT

    # Numero de Knudsen
    lam = m0 / (np.sqrt(2) * np.pi * sigma ** 2 * rho * Na)
    Kn = lam / lc

    # Parametro de regime de escoamento
    if Kn >= 10:
        d = 1
    elif Kn <= 0.01:
        d = 2
    else:
        d = 3

    # Re - Numero de Reynolds [adm]
    Re = rho * v * lc / mu

    return T, Tm, p, rho, ainf, M, mu, Pr, Kn, d, Re, R
