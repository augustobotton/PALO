import numpy as np
import json
import os


class ModeloAtmosferico:
    def __init__(self, json_file_path):
        # Verifica se o arquivo existe no caminho especificado
        if not os.path.exists(json_file_path):
            raise FileNotFoundError(f"O arquivo {json_file_path} não foi encontrado.")

        # Carrega os dados do arquivo JSON
        with open(json_file_path, 'r') as file:
            data = json.load(file)

        # Converte os dados do JSON em arrays numpy
        self.hi = np.array(data['hi']) * 1e3  # km para m
        self.Ti = np.array(data['Ti'])
        self.Ri = np.array(data['Ri'])
        self.a = np.array(data['a']) * 1e-3  # °C/km para °C/m
        self.Mi = np.array(data['Mi'])

        # Carrega as constantes do JSON
        self.g0 = data['constantes']['aceleracao_gravidade']  # Aceleração da gravidade ao nível do mar (m/s²)
        self.Na = data['constantes']['numero_avogadro']  # Número de Avogadro
        self.sigma = data['constantes']['diametro_colisao']  # Diâmetro de colisão para o ar (m)
        self.m0 = data['constantes']['massa_molar_ar']  # Massa molar do ar ao nível do mar (kg/mol)
        self.P0 = data['constantes']['pressao_nivel_mar']  # Pressão ao nível do mar (N/m²)
        self.Requat = data['constantes']['raio_medio_terra']  # Raio médio da Terra (m)
        self.gamma = data['constantes'][
            'razao_calores_especificos']  # Razão de calores específicos ao nível do mar
        self.beta = 2 / self.Requat  # Número beta associado à distância radial média do nível do mar

    def calcula(self, altitude_geometrica, v, lc, dT):
        if altitude_geometrica < 0:
            i = 0
            altitude_geometrica = 0
        elif altitude_geometrica > 2000e3:
            i = 20
            altitude_geometrica = 2000e3
        else:
            for i in range(21):
                if i == 20 or (altitude_geometrica >= self.hi[i] and altitude_geometrica < self.hi[i + 1]):
                    break

        # Calcula a pressão inicial da camada i
        Pi = self.P0
        px = self.P0
        for j in range(1, i + 1):
            if self.a[j - 1] != 0:
                A = 1 + (self.a[j - 1] * (self.hi[j] - self.hi[j - 1])) / self.Ti[j - 1]
                B = -(self.g0 / (self.Ri[j] * self.a[j - 1])) * (
                        1 + self.beta * ((self.Ti[j - 1] / self.a[j - 1]) - self.hi[j - 1]))
                C = (self.g0 * self.beta / (self.Ri[j] * self.a[j - 1])) * (self.hi[j] - self.hi[j - 1])
                Pi = px * (A ** B) * np.exp(C)
                px = Pi
            else:
                Pi = px * np.exp(
                    -(self.g0 / (self.Ri[j] * self.Ti[j - 1])) * (self.hi[j] - self.hi[j - 1]) * (
                            1 - (self.beta / 2) * (self.hi[j] - self.hi[j - 1])))
                px = Pi

        temperatura_escala_molecular = self.Ti[i] + self.a[i] * (altitude_geometrica - self.hi[i]) + dT
        if i >= 20:
            constante_de_gas_ideal = self.Ri[i]
            M = self.Mi[i]
        else:
            constante_de_gas_ideal = self.Ri[i] + (
                        (self.Ri[i + 1] - self.Ri[i]) / (self.hi[i + 1] - self.hi[i])) * (
                                             altitude_geometrica - self.hi[i])
            M = self.Mi[i] + ((self.Mi[i + 1] - self.Mi[i]) / (self.hi[i + 1] - self.hi[i])) * (
                    altitude_geometrica - self.hi[i])

        temperatura_cinetica = (M / self.m0) * temperatura_escala_molecular
        if self.a[i] != 0:
            A = 1 + (self.a[i] * (altitude_geometrica - self.hi[i])) / self.Ti[i]
            B = -(self.g0 / (constante_de_gas_ideal * self.a[i])) * (
                    1 + self.beta * ((self.Ti[i] / self.a[i]) - self.hi[i]))
            C = (self.g0 * self.beta / (constante_de_gas_ideal * self.a[i])) * (
                        altitude_geometrica - self.hi[i])
            pressao = Pi * (A ** B) * np.exp(C)
        else:
            pressao = Pi * np.exp(
                -(self.g0 / (constante_de_gas_ideal * self.Ti[i])) * (altitude_geometrica - self.hi[i]) * (
                        1 - (self.beta / 2) * (altitude_geometrica - self.hi[i])))

        densidade_do_ar = pressao / (constante_de_gas_ideal * temperatura_escala_molecular)
        velocidade_do_som = np.sqrt(self.gamma * constante_de_gas_ideal * temperatura_escala_molecular)
        numero_de_mach = v / velocidade_do_som
        viscosidade_dinamica = 1.458e-6 * (temperatura_escala_molecular ** (3 / 2)) / (
                temperatura_escala_molecular + 110.4)
        cp = constante_de_gas_ideal * self.gamma / (self.gamma - 1)
        kT = (2.64638e-3 * (temperatura_escala_molecular ** (3 / 2))) / (
                temperatura_escala_molecular + 245.4 * (10 ** (-12 / temperatura_escala_molecular)))
        numero_de_prandtl = viscosidade_dinamica * cp / kT
        lam = self.m0 / (np.sqrt(2) * np.pi * self.sigma ** 2 * densidade_do_ar * self.Na)
        numero_de_knudsen = lam / lc
        if numero_de_knudsen >= 10:
            parametro_regime_de_escoamento = 1
        elif numero_de_knudsen <= 0.01:
            parametro_regime_de_escoamento = 2
        else:
            parametro_regime_de_escoamento = 3
        reynolds = densidade_do_ar * v * lc / viscosidade_dinamica

        return temperatura_cinetica, temperatura_escala_molecular, pressao, densidade_do_ar, velocidade_do_som, numero_de_mach, viscosidade_dinamica, numero_de_prandtl, numero_de_knudsen, parametro_regime_de_escoamento, reynolds,constante_de_gas_ideal


