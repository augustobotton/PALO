import numpy as np

#TODO trocar das rotações para os angulos de acordo com as referencias

# Matrizes de rotação elementares
def rotX(alf):
	"""
    Gera a matriz de rotação ao redor do eixo X por um ângulo alf.

    Parâmetros:
    alf (float): Ângulo de rotação em radianos.

    Retorna:
    numpy.ndarray: Matriz de rotação 3x3.
    """
	C1 = np.array([
		[1, 0, 0],
		[0, np.cos(alf), np.sin(alf)],
		[0, -np.sin(alf), np.cos(alf)]
	])
	return C1


def rotY(alf):
	"""
    Gera a matriz de rotação ao redor do eixo Y por um ângulo alf.

    Parâmetros:
    alf (float): Ângulo de rotação em radianos.

    Retorna:
    numpy.ndarray: Matriz de rotação 3x3.
    """
	C2 = np.array([
		[np.cos(alf), 0, -np.sin(alf)],
		[0, 1, 0],
		[np.sin(alf), 0, np.cos(alf)]
	])
	return C2


def rotZ(alf):
	"""
    Gera a matriz de rotação ao redor do eixo Z por um ângulo alf.

    Parâmetros:
    alf (float): Ângulo de rotação em radianos.

    Retorna:
    numpy.ndarray: Matriz de rotação 3x3.
    """
	C3 = np.array([
		[np.cos(alf), np.sin(alf), 0],
		[-np.sin(alf), np.cos(alf), 0],
		[0, 0, 1]
	])
	return C3


def matInercPerif(OMEGA, inc, omega):
	"""
    Calcula a matriz de rotação do sistema de referência inercial para o sistema perifocal.

    Parâmetros:
    OMEGA (float): Longitude do nó ascendente em radianos.
    inc (float): Inclinação da órbita em radianos.
    omega (float): Argumento do perigeu em radianos.

    Retorna:
    numpy.ndarray: Matriz de rotação 3x3 do sistema inercial para o perifocal.
    """
	Cip = rotZ(omega) @ rotX(inc) @ rotZ(OMEGA)
	return Cip
