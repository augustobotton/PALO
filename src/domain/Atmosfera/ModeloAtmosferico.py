import numpy as np


class ModeloAtmosferico:
    def __init__(self):

        self.hi = 1e3 * np.array([0, 11.0191, 20.0631, 32.1619, 47.3501, 51.4125, 71.8020, 86, 100, 110, 120, 150,
                                  160, 170, 190, 230, 300, 400, 500, 600, 700])  # km
        # m
        self.Ti = np.array([288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65, 186.946, 210.65, 260.65, 360.65,
                            960.65, 1110.65, 1210.65, 1350.65, 1550.65, 1830.65, 2160.65, 2420.65, 2590.65, 2700.65])
        self.Ri = np.array([287.0, 287.0, 287.0, 287.0, 287.0, 287.0, 287.02, 287.02, 287.84, 291.06, 308.79, 311.80,
                            313.69, 321.57, 336.68, 366.84, 416.88, 463.36, 493.63, 514.08, 514.08])
        self.a = 1e-3 * np.array([-6.5, 0.0, 1.0, 2.8, 0.0, -2.8, -2.0, 1.693, 5.0, 10.0, 20.0, 15.0, 10.0, 7.0, 5.0,
                                  4.0, 3.3, 2.6, 1.7, 1.1, 0.0])  # °/km

        self.Mi = np.array([28.9644, 28.9644, 28.9644, 28.9644, 28.9644, 28.9644, 28.9644, 28.9644,
                            28.88, 28.56, 28.07, 26.92, 26.66, 26.4, 25.85, 24.7, 22.66, 19.94, 17.94, 16.84, 16.17])

        # Constantes
        self.M0 = self.Mi[0]  # Peso molecular ao nivel do mar
        self.g0 = 9.80665  # Valor ao nivel do mar da aceleracao da gravidade (m/s**2)
        self.Na = 6.0220978e23  # Numero de Avogadro
        self.sigma = 3.65e-10  # Diametro de colisao para o ar (m)
        self.m0 = 28.964e-3  # Massa molar do ar ao nivel do mar (kg/Mol)
        self.P0 = 1.01325e5  # Pressao padrao ao nivel do mar (N/m**2)
        self.Requat = 6378.14e3  # Raio medio da Terra (m)
        self.gamma = 1.405  # Razao de calores especificos ao nivel do mar
        # Constantes calculadas
        # Numero beta associado a distancia radial media do nivel do mar
        self.beta = 2 / self.Requat
        self.temperatura_cinetica = 0.
        self.temperatura_escala_molecular = 0.
        self.pressao = 0.
        self.densidade_do_ar = 0.
        self.velocidade_do_som = 0.
        self.numero_de_mach = 0.
        self.viscosidade_dinamica = 0.
        self.numero_de_prandtl = 0.
        self.numero_de_knudsen = 0.
        self.parametro_regime_de_escoamento = 0.
        self.reynolds = 0.
        self.constante_de_gas_ideal = 0.

    def calcula(self, altitude_geometrica, v, lc, dT):

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
                elif ((altitude_geometrica >= self.hi[i]) and (altitude_geometrica < self.hi[i + 1])):
                    break

        # Realiza os calculos
        # Pressao no inicio da camada i
        Pi = self.P0  # Pressao inicial da camada - inicializa com o valor ao nivel do mar
        px = self.P0  # Variavel auxiliar para guardar o valor de pressao inicial na camada anterior
        for j in range(1, i + 1):
            if self.a[j - 1] != 0:  # Verifica se a camada nao eh isotermica
                # Calcula a pressao inicial da camada j a partir do modelo da
                # camada j-1
                A = 1 + (self.a[j - 1] * (self.hi[j] - self.hi[j - 1])) / self.Ti[j - 1]
                B = -(self.g0 / (self.Ri[j] * self.a[j - 1])) * (
                        1 + self.beta * ((self.Ti[j - 1] / self.a[j - 1]) - self.hi[j - 1]))
                C = (self.g0 * self.beta / (self.Ri[j] * self.a[j - 1])) * (self.hi[j] - self.hi[j - 1])
                Pi = px * (A ** B) * np.exp(C)
                px = Pi  # O valor atual sera o valor anterior na proxima iteracao
            else:
                # Calcula a pressao inicial da camada i pelo modelo isotermico
                Pi = px * np.exp(
                    -(self.g0 / (self.Ri[j] * self.Ti[j - 1])) * (self.hi[j] - self.hi[j - 1]) * (
                            1 - (self.beta / 2) * (self.hi[j] - self.hi[j - 1])))
                px = Pi  # O valor atual sera o valor anterior na proxima iteracao
        # Temperatura padrao (molecular) - usada nos calculos internos
        temperatura_escala_molecular = self.Ti[i] + self.a[i] * (altitude_geometrica - self.hi[i])
        temperatura_escala_molecular = temperatura_escala_molecular + dT

        # Interpola o valor de R e M

        if i >= 20:
            constante_de_gas_ideal = self.Ri[i]
            M = self.Mi[i]
        else:
            constante_de_gas_ideal = self.Ri[i] + ((self.Ri[i + 1] - self.Ri[i]) / (self.hi[i + 1] - self.hi[i])) * (
                    altitude_geometrica - self.hi[i])
            M = self.Mi[i] + ((self.Mi[i + 1] - self.Mi[i]) / (self.hi[i + 1] - self.hi[i])) * (
                    altitude_geometrica - self.hi[i])
            # Temperatura padrao (cinetica) - saida da funcao
        temperatura_cinetica = (M / self.M0) * temperatura_escala_molecular
        # Pressao
        if self.a[i] != 0:  # Verifica se a camada nao eh isotermica
            # Calcula a pressao
            A = 1 + (self.a[i] * (altitude_geometrica - self.hi[i])) / self.Ti[i]
            B = -(self.g0 / (constante_de_gas_ideal * self.a[i])) * (
                    1 + self.beta * ((self.Ti[i] / self.a[i]) - self.hi[i]))
            C = (self.g0 * self.beta / (constante_de_gas_ideal * self.a[i])) * (altitude_geometrica - self.hi[i])
            pressao = Pi * (A ** B) * np.exp(C)
        else:
            # Calcula a pressao pelo modelo isotermico
            pressao = Pi * np.exp(
                -(self.g0 / (constante_de_gas_ideal * self.Ti[i])) * (altitude_geometrica - self.hi[i]) * (
                        1 - (self.beta / 2) * (altitude_geometrica - self.hi[i])))

        # Calcula a densidade
        densidade_do_ar = pressao / (constante_de_gas_ideal * temperatura_escala_molecular)
        # Velocidade do som
        velocidade_do_som = np.sqrt(self.gamma * constante_de_gas_ideal * temperatura_escala_molecular)
        # Numero de numero_de_mach
        numero_de_mach = v / velocidade_do_som
        # Coeficiente de viscosidade dinamica
        viscosidade_dinamica = 1.458e-6 * (temperatura_escala_molecular ** (3 / 2)) / (
                temperatura_escala_molecular + 110.4)
        # Numero de Prandtl
        cp = constante_de_gas_ideal * self.gamma / (self.gamma - 1)
        kT = (2.64638e-3 * (temperatura_escala_molecular ** (3 / 2))) / (
                temperatura_escala_molecular + 245.4 * (10 ** (-12 / temperatura_escala_molecular)))
        numero_de_prandtl = viscosidade_dinamica * cp / kT
        # Numero de Knudsen
        lam = self.m0 / (np.sqrt(2) * np.pi * self.sigma ** 2 * densidade_do_ar * self.Na)
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
