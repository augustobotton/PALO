import numpy as np
from scipy.optimize import fsolve

from src.domain.utilidades_mecanica_orbital.propagacao import sistReferencia
from src.domain.utilidades_mecanica_orbital.Orbitas.ConstrutorOrbita import ConstrutorOrbita
from src.domain.utilidades_mecanica_orbital.Orbitas.ModeloOrbita import Orbita


#
# Módulo com funções de propagação de órbita kepleriana


# Órbita elíptica: cálculo da propagação temporal da anomalia excêntrica
#


def resolveEqKepler(t, orbita: Orbita):
	# Órbita eliptica: calcula a anomalia excentrica para dado tempo
	# Entrada
	# t: [s] tempo
	# mu: [m^3/s^2] parâmetro gravitacional do primário
	# orb: vetor de elementos orbitais clássicos
	# orb(0): a [m] - semi eixo maior
	# orb(1): e - excentricidade
	# orb(2): tau [s] - tempo de periastro, medido com respeito a t=0
	# Saida
	# E: [rad] anomalia excentrica no instante t
	# Elementos orbitais
	a = orbita.semi_eixo_maior
	e = orbita.excentricidade  # Excentricidade
	mu = orbita.mu
	tau = orbita.tempo_de_periastro  # tempo de periastro

	# Movimento medio
	n = np.sqrt(mu / a ** 3)

	# Anomalia media
	M = n * (t - tau)

	# Passagem de parâmetros para a função objetivo
	dados = (e, M)

	# Chute inicial
	E0 = M  # Seria o resultado da orbita circular

	E = fsolve(Kepler, E0, args=dados)
	return E[0]


# fsolve é bandida, temos que indexar para obter a solução em float, não ela
# retorna uma tupla, mesmo que a solução seja só um valor.
#
# Função para propagação de órbita elíptica no referencial inercial


def propagaEliptica(t, orbita):
	# Entradas:
	# t: [s] tempo. A sua referência é o tempo de periastro fornecido.
	# mu: [m^3/s^2] parâmetro gravitacional do primário.
	# orb: vetor 6x1 de elementos orbitais clássicos
	# orb[0]: a [m] - semi eixo maior
	# orb[1]: e - excentricidade
	# orb[2]: tau [s] - tempo de periastro, medido em relação ao tempo t=0
	# orb[3]: OMEGA [rad] - longitude celeste do nodo ascendente do referencial
	# perifocal com respeito ao inercial.
	# orb[4]: i [rad] - inclinação da órbita com respeito ao plano XY do referencialinercial
	# orb[5]: omega [rad] - argumento de periastro
	# Saídas:
	# Variáveis calculadas no instante t fornecido
	# theta: [rad] - anomalia verdadeira
	# Ri: [m] - Vetor 3x1, posição no referencial inercial. Coordenadas retangulares
	# Vi: [m] - Vetor 3x1, velocidade no referencial inercial. Coordenadas retangulares

	# Elementos orbitais
	a = orbita.semi_eixo_maior
	e = orbita.excentricidade
	OMEGA = orbita.raan
	inc = orbita.inclinacao
	omega = orbita.arg_periastro
	mu = orbita.mu

	# Anomalia verdadeira em t=0
	E0 = resolveEqKepler(0, orbita)

	theta0 = 2 * np.arctan(np.sqrt((1 + e) / (1 - e)) * np.tan(E0 / 2))

	# Parametro
	p = orbita.calcular_parametro_orbital()
	# Período
	n = np.sqrt(mu / a ** 3);  # movimento médio da orbita elíptica
	f = n / (2 * np.pi);  # frequência
	P = 1 / f  # Período
	# Vetores posição e velocidade inicial, no referencial perifocal, escritos
	# em coordenadas retangulares
	h = np.sqrt(p * mu)  # [m^2/s] Quantidade de movimento angular especifica
	r0 = p / (1 + e * np.cos(theta0))
	R0 = np.array([r0 * np.cos(theta0), r0 * np.sin(theta0), 0])
	V0 = np.array([-(mu / h) * np.sin(theta0), (mu / h) * (e + np.cos(theta0)), 0])
	# Resolve a equacao de Kepler, determinando a anomalia excêntrica
	E = resolveEqKepler(t, orbita)
	# Determina a anomalia verdadeira
	theta = 2 * np.arctan(np.sqrt((1 + e) / (1 - e)) * np.tan(E / 2))
	if theta < 0:
		theta = theta + 2 * np.pi
	# Determina a matriz de transição de estado a partir de theta
	PHI = matrizTransicaoEstado(theta, theta0, orbita)
	# Determina posição e velocidade perifocal em função da anomalia
	# verdadeira pela matriz de transição de estado
	R = PHI[0, 0] * R0 + PHI[0, 1] * V0
	V = PHI[1, 0] * R0 + PHI[1, 1] * V0
	# Matriz de transformação de coordenadas do referencial inercial para o perifocal
	Cip = sistReferencia.matInercPerif(OMEGA, inc, omega)
	# Do perifocal para o inercial
	Cpi = np.transpose(Cip)
	# Posição e velocidade no referencial inercial
	Ri = Cpi @ R;
	Vi = Cpi @ V
	return theta, Ri, Vi


def Kepler(E, *dados):
	# Funcao objetivo: Equação de Kepler
	# Entradas:
	# E: [rad] anomalia excêntrica - incógnita
	# *dados: parâmetros da função
	# Saída:
	# y: quando igual a zero, a equacao de Kepler esta resolvida
	#
	# Recebimento de parametros
	e, M = dados

	y = E - e * np.sin(E) - M
	return y


#


def resolveEqBarker(t: float, orbita: Orbita) -> float:
	"""
  Calcula a anomalia verdadeira - usando a solução analitica da equação de Barker - para uma órbita parabólica
  em um dado momento.

  Parâmetros:
  t (float): tempo desde o periastro, em segundos.
  orbita (Orbita): objeto da classe Orbita representando a órbita parabólica.
  mu (float): parâmetro gravitacional padrão (GM) do corpo central.

  Retorna:
  float: anomalia verdadeira, em radianos.
  """
	# Cálculo do parâmetro orbital
	parametro_orbital = orbita.parametro

	# Anomalia média parabólica
	mp = np.sqrt(orbita.mu / parametro_orbital ** 3) * (t - orbita.tempo_de_periastro)

	# Solução analítica para tangente de meio theta
	termo_raiz = np.sqrt(1 + 9 * mp ** 2)
	tan_meio_theta = (3 * mp + termo_raiz) ** (1 / 3) - (3 * mp + termo_raiz) ** (-1 / 3)

	# Cálculo da anomalia verdadeira
	theta = 2 * np.arctan(tan_meio_theta)
	return theta


def KeplerHiperbolica(H, *dados):
	# Função objetivo: equação de Kepler hiperbólica
	# Entrada
	# H: [rad] anomalia hiperbólica
	# *dados: parâmetros da função
	# Saida
	# y: quando igual a zero, a equacao de Kepler hiperbólica está resolvida
	#
	# Recebimento de parametros
	e, Mh = dados
	# Quando y=0, a equacao de Kepler hiperbólica está resolvida
	y = e * np.sinh(H) - H - Mh
	return y


def resolveEqKeplerHiperbolica(tempo: float, orbita: Orbita):
	# Órbita hiperbólica: calcula a anomalia hiperbólica para dado tempo
	# Entrada
	# tempo: [s] tempo
	# Saida

	# H: [rad] anomalia hiperbólica

	tau = orbita.tempo_de_periastro  # tempo de periastro
	p = orbita.calcular_parametro_orbital()

	e = orbita.excentricidade  # Excentricidade da orbita
	# Anomalia media hiperbolica
	Mh = (e ** 2 - 1) ** (3 / 2) * np.sqrt(orbita.mu / p ** 3) * (tempo - tau)
	# Passagem de parâmetros para a função objetivo
	dados = (e, Mh)
	# Chute inicial
	H0 = Mh

	H = fsolve(KeplerHiperbolica, H0, args=dados)
	return H


def matrizTransicaoEstado(theta, theta0, orbita):
	# Funcao para determinar a matriz de transicao de estado cujos elementos
	# sao os coeficientes de Lagrange
	# Entradas:
	# theta [rad]: anomalia verdadeira
	# mu: [m^3/s^2] parâmetro gravitacional do primário.
	# Saída:
	# PHI [2x2]: matriz de transicao de estado

	# Elementos orbitais
	e = orbita.excentricidade
	mu = orbita.mu
	# Calculo das constantes associadas a orbita
	p = orbita.calcular_parametro_orbital()
	h = np.sqrt(p * mu)  # momento angular especifico da orbita
	r0 = p / (1 + e * np.cos(theta0))  # distancia radial inicial

	# Calculo dos coeficientes de Lagrange
	r = p / (1 + e * np.cos(theta))  # Distancia radial para o theta dado
	f = 1 + (r / p) * (np.cos(theta - theta0) - 1)
	g = (r * r0 / h) * np.sin(theta - theta0)
	fp = -(h / p ** 2) * (np.sin(theta - theta0) + e * (np.sin(theta) - np.sin(theta0)))
	gp = 1 + (r0 / p) * (np.cos(theta - theta0) - 1)

	# Matriz de transicao de estado
	return np.array([[f, g], [fp, gp]])






