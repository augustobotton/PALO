import numpy as np


def atm_padrao(altitude_geometrica, v, lc, dT):
    # Funcao para calculo da atmosfera padrao para a altitude geometrica de
    # Válido de 0 a 200 km
    # É utilizado o modelo apresentado na referencia
    # TEWARI, A. Atmospheric and Space Flight Dynamics:
    # Modelling and simulation with MATLAB and Simulink. Boston: Birkhauser, 2007.
    # Capitulo 9 - Planetary atmosphere
    # O modelo possui 21 camadas, sendo uma combinacao dos modelos de atmosfera
    # padrao norte americanos de 1962 e 1976.
    # Todas as camadas consideradas possuem variacao linear de temperatura ou
    # sao isotermicas. Na faixa de altitude de zero ate 86km, eh utilizado o
    # perfil de temperaturas da atmosfera padrao norte americana do ano de
    # 1976. Para a faixa de altitude de 86km a 2000km, eh utilizado o modelo de 1962.
    # Apesar de a atmosfera nao apresentar equilibrio termico e quimico acima
    # de 86km, sendo os perfis de temperatura, nessas regioes, nao lineares,
    # esta funcao adota do modelo linear de temperaturas de 1962 como uma aproximacao.
    # Entradas:
    # h - altitude geometrica [m]
    # v - velocidade do veiculo em relacao ao escoamento nao perturbado, em
    # metros por segundo [m/s]
    # lc - comprimento caracteristico do veiculo, em metros [m]
    # Saidas
    # T - Temperatura cinetica, em Kelvin [K]
    # Tm - Temperatura de escala molecular, em Kelvin [K]
    # p - Pressao, em Pascal [N/m**2]
    # rho - Densidade [kg/m**2]
    # ainf - Velocidade do som [m/s]
    # M - Numero de numero_de_mach [adm]
    # mu - Coeficiente de viscosidade dinamica [kg/m*s]
    # Pr - Numero de Prandtl [adm]
    # Kn - Numero de Knudsen [adm]
    # d - Parametro de regime de escoamento [adm]
    # reynolds - Numero de Reynolds [adm]

    # Entrada de dados
    # Vetores com os seguintes dados, cada linha eh uma camada do modelo de
    # atmosfera padrao: altitude no inicio da camada (m), temperatura no inicio
    # da camada (K), constante de gas ideal do ar na camada (J/kg.K), taxa de
    # lapso termico (K/m), peso molecular M

    hi = 1e3 * np.array([0, 11.0191, 20.0631, 32.1619, 47.3501, 51.4125, 71.8020, 86, 100, 110, 120, 150,
                         160, 170, 190, 230, 300, 400, 500, 600, 700])  # km
    # m
    Ti = np.array([288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65, 186.946, 210.65, 260.65, 360.65,
                   960.65, 1110.65, 1210.65, 1350.65, 1550.65, 1830.65, 2160.65, 2420.65, 2590.65, 2700.65])
    Ri = np.array([287.0, 287.0, 287.0, 287.0, 287.0, 287.0, 287.02, 287.02, 287.84, 291.06, 308.79, 311.80,
                   313.69, 321.57, 336.68, 366.84, 416.88, 463.36, 493.63, 514.08, 514.08])
    a = 1e-3 * np.array([-6.5, 0.0, 1.0, 2.8, 0.0, -2.8, -2.0, 1.693, 5.0, 10.0, 20.0, 15.0, 10.0, 7.0, 5.0,
                         4.0, 3.3, 2.6, 1.7, 1.1, 0.0])  # °/km

    Mi = np.array([28.9644, 28.9644, 28.9644, 28.9644, 28.9644, 28.9644, 28.9644, 28.9644,
                   28.88, 28.56, 28.07, 26.92, 26.66, 26.4, 25.85, 24.7, 22.66, 19.94, 17.94, 16.84, 16.17])

    # Constantes
    M0 = Mi[0]  # Peso molecular ao nivel do mar
    g0 = 9.80665  # Valor ao nivel do mar da aceleracao da gravidade (m/s**2)
    Na = 6.0220978e23  # Numero de Avogadro
    sigma = 3.65e-10  # Diametro de colisao para o ar (m)
    m0 = 28.964e-3  # Massa molar do ar ao nivel do mar (kg/Mol)
    P0 = 1.01325e5  # Pressao padrao ao nivel do mar (N/m**2)
    Requat = 6378.14e3  # Raio medio da Terra (m)
    gamma = 1.405  # Razao de calores especificos ao nivel do mar
    # Constantes calculadas
    # Numero beta associado a distancia radial media do nivel do mar
    beta = 2 / Requat

    # Identifica a camada a qual a altitude pertence

    # Contador
    i = 0
    if altitude_geometrica < 0:
        # Altitude negativa. Os resultados apresentados dizem respeito a h=0.
        i = 0;
        altitude_geometrica = 0
    elif altitude_geometrica > 2000e3:
        # A altitude fornecida esta acima do limite superior de 2.000 km.
        # Os resultados apresentados dizem respeito a h=2.000 km.
        i = 20;
        altitude_geometrica = 2000e3
    else:
        for i in range(21):
            if (i == 20):
                break
            elif ((altitude_geometrica >= hi[i]) and (altitude_geometrica < hi[i + 1])):
                break

    # Realiza os calculos
    # Pressao no inicio da camada i
    Pi = P0  # Pressao inicial da camada - inicializa com o valor ao nivel do mar
    px = P0  # Variavel auxiliar para guardar o valor de pressao inicial na camada anterior
    for j in range(1, i + 1):
        if a[j - 1] != 0:  # Verifica se a camada nao eh isotermica
            # Calcula a pressao inicial da camada j a partir do modelo da
            # camada j-1
            A = 1 + (a[j - 1] * (hi[j] - hi[j - 1])) / Ti[j - 1]
            B = -(g0 / (Ri[j] * a[j - 1])) * (1 + beta * ((Ti[j - 1] / a[j - 1]) - hi[j - 1]))
            C = (g0 * beta / (Ri[j] * a[j - 1])) * (hi[j] - hi[j - 1])
            Pi = px * (A ** B) * np.exp(C)
            px = Pi  # O valor atual sera o valor anterior na proxima iteracao
        else:
            # Calcula a pressao inicial da camada i pelo modelo isotermico
            Pi = px * np.exp(-(g0 / (Ri[j] * Ti[j - 1])) * (hi[j] - hi[j - 1]) * (1 - (beta / 2) * (hi[j] - hi[j - 1])))
            px = Pi  # O valor atual sera o valor anterior na proxima iteracao
    # Temperatura padrao (molecular) - usada nos calculos internos
    temperatura_escala_molecular = Ti[i] + a[i] * (altitude_geometrica - hi[i])
    temperatura_escala_molecular = temperatura_escala_molecular + dT

    # Interpola o valor de R e M

    if i >= 20:
        constante_de_gas_ideal = Ri[i]
        M = Mi[i]
    else:
        constante_de_gas_ideal = Ri[i] + ((Ri[i + 1] - Ri[i]) / (hi[i + 1] - hi[i])) * (altitude_geometrica - hi[i])
        M = Mi[i] + ((Mi[i + 1] - Mi[i]) / (hi[i + 1] - hi[i])) * (altitude_geometrica - hi[i])
        # Temperatura padrao (cinetica) - saida da funcao
    temperatura_cinetica = (M / M0) * temperatura_escala_molecular
    # Pressao
    if a[i] != 0:  # Verifica se a camada nao eh isotermica
        # Calcula a pressao
        A = 1 + (a[i] * (altitude_geometrica - hi[i])) / Ti[i]
        B = -(g0 / (constante_de_gas_ideal * a[i])) * (1 + beta * ((Ti[i] / a[i]) - hi[i]))
        C = (g0 * beta / (constante_de_gas_ideal * a[i])) * (altitude_geometrica - hi[i])
        pressao = Pi * (A ** B) * np.exp(C)
    else:
        # Calcula a pressao pelo modelo isotermico
        pressao = Pi * np.exp(-(g0 / (constante_de_gas_ideal * Ti[i])) * (altitude_geometrica - hi[i]) * (
                    1 - (beta / 2) * (altitude_geometrica - hi[i])))

    # Calcula a densidade
    densidade_do_ar = pressao / (constante_de_gas_ideal * temperatura_escala_molecular)
    # Velocidade do som
    velocidade_do_som = np.sqrt(gamma * constante_de_gas_ideal * temperatura_escala_molecular)
    # Numero de numero_de_mach
    numero_de_mach = v / velocidade_do_som
    # Coeficiente de viscosidade dinamica
    viscosidade_dinamica = 1.458e-6 * (temperatura_escala_molecular ** (3 / 2)) / (temperatura_escala_molecular + 110.4)
    # Numero de Prandtl
    cp = constante_de_gas_ideal * gamma / (gamma - 1)
    kT = (2.64638e-3 * (temperatura_escala_molecular ** (3 / 2))) / (
                temperatura_escala_molecular + 245.4 * (10 ** (-12 / temperatura_escala_molecular)))
    numero_de_prandtl = viscosidade_dinamica * cp / kT
    # Numero de Knudsen
    lam = m0 / (np.sqrt(2) * np.pi * sigma ** 2 * densidade_do_ar * Na)
    numero_de_knudsen = lam / lc
    # Parametro de regime de escoamento
    if numero_de_knudsen >= 10:
        parametro_regime_de_escoamento = 1
    elif numero_de_knudsen <= 0.01:
        parametro_regime_de_escoamento = 2
    else:
        parametro_regime_de_escoamento = 3
    # Re - Numero de Reynolds [adm]
    reynolds = densidade_do_ar * v * lc / viscosidade_dinamica
    return temperatura_cinetica, temperatura_escala_molecular, pressao, densidade_do_ar, velocidade_do_som, numero_de_mach, viscosidade_dinamica, numero_de_prandtl, numero_de_knudsen, parametro_regime_de_escoamento, reynolds, constante_de_gas_ideal
