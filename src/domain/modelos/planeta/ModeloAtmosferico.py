import json
import os

import numpy as np


class ModeloAtmosferico:
    """
    Modelo atmosférico baseado em dados de um arquivo JSON que define camadas atmosféricas e constantes físicas.

    Attributes:
        hi (np.ndarray): Alturas das camadas atmosféricas em metros.
        Ti (np.ndarray): Temperaturas de base para cada camada em Kelvin.
        Ri (np.ndarray): Constante de gás ideal para cada camada em J/(kg·K).
        a (np.ndarray): Gradiente térmico em °C por metro.
        Mi (np.ndarray): Massa molar média para cada camada em kg/mol.
        g0 (float): Aceleração da gravidade ao nível do mar em m/s².
        Na (float): Número de Avogadro.
        sigma (float): Diâmetro de colisão molecular em metros.
        m0 (float): Massa molar do ar ao nível do mar em kg/mol.
        P0 (float): Pressão atmosférica ao nível do mar em N/m².
        Requat (float): Raio médio da Terra em metros.
        gamma (float): Razão de calores específicos.
        beta (float): Coeficiente relacionado ao raio médio da Terra.
    """

    def __init__(self, json_file_path):
        """
        Inicializa o ModeloAtmosferico carregando dados de um arquivo JSON.

        Args:
            json_file_path (str): Caminho do arquivo JSON com dados atmosféricos.

        Raises:
            FileNotFoundError: Se o arquivo JSON não for encontrado.
        """
        if not os.path.exists(json_file_path):
            raise FileNotFoundError(f"O arquivo {json_file_path} não foi encontrado.")

        with open(json_file_path, 'r') as file:
            data = json.load(file)

        self.hi = np.array(data['hi']) * 1e3  # Convertendo km para m
        self.Ti = np.array(data['Ti'])
        self.Ri = np.array(data['Ri'])
        self.a = np.array(data['a']) * 1e-3  # Convertendo °C/km para °C/m
        self.Mi = np.array(data['Mi'])

        self.g0 = data['constantes']['aceleracao_gravidade']
        self.Na = 6.0220978e23
        self.sigma = 3.65e-10
        self.m0 = self.Mi[0]
        self.P0 = data['constantes']['pressao_nivel_mar']
        self.Requat = data['constantes']['raio_medio_terra']
        self.gamma = data['constantes']['razao_calores_especificos']
        self.beta = 2 / self.Requat

    def calcula(self, altitude_geometrica, v, lc, dT):

        """
                Calcula e retorna várias propriedades atmosféricas em uma dada altitude.

                Args:
                    altitude_geometrica (float): Altitude geométrica em metros.
                    v (float): Velocidade em m/s.
                    lc (float): Comprimento característico em metros.
                    dT (float): Variação de temperatura local em Kelvin.

                Returns:
                    tuple: Contém temperatura cinética, temperatura escala molecular, pressão, densidade do ar,
                           velocidade do som, número de Mach, viscosidade dinâmica, número de Prandtl,
                           número de Knudsen, parâmetro regime de escoamento, número de Reynolds e constante de gás ideal.
                """

        if altitude_geometrica < 0:
            i = 0
            altitude_geometrica = 0
        elif altitude_geometrica > 2000e3:
            i = 20
            altitude_geometrica = 2000e3
        else:
            for i in range(21):
                if (i == 20):
                    break
                elif ((altitude_geometrica >= self.hi[i]) and (altitude_geometrica < self.hi[i + 1])):
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

        temperatura_escala_molecular = self.Ti[i] + self.a[i] * (altitude_geometrica - self.hi[i]) 
        temperatura_escala_molecular = temperatura_escala_molecular +dT
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
        numero_de_knudsen = lam / lc/1000
        if numero_de_knudsen >= 10:
            parametro_regime_de_escoamento = 1
        elif numero_de_knudsen <= 0.01:
            parametro_regime_de_escoamento = 2
        else:
            parametro_regime_de_escoamento = 3
        reynolds = densidade_do_ar * v * lc / viscosidade_dinamica

        return temperatura_cinetica, temperatura_escala_molecular, pressao, densidade_do_ar, velocidade_do_som, numero_de_mach, viscosidade_dinamica, numero_de_prandtl, numero_de_knudsen, parametro_regime_de_escoamento, reynolds,constante_de_gas_ideal


