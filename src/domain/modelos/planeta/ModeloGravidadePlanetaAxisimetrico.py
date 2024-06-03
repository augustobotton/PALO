import numpy as np


def calcular_gravidade_axisimetric(distancia_radial, colatitude, planeta):
	"""
    Calcula a gravidade de um planeta axissimétrico considerando as constantes de Jeffery (J2, J3, J4).

    Parâmetros:
    distancia_radial (float): Distância radial ao centro de massa do planeta (m).
    colatitude (float): Colatitude (rad).
    planeta (object): Objeto contendo os parâmetros raio_equatorial, gravidade_padrao_nivel_do_mar,
                      constante_gravitacional, J2, J3, J4.

    Retorna:
    tuple: Componente radial da gravidade (m/s^2) e componente colatitudinal (m/s^2).
    """
	# Calculando a componente radial da gravidade
	componente_radial = (
			(1 / distancia_radial ** 6)
			* planeta.gravidade_padrao_nivel_do_mar
			* planeta.constante_gravitacional
			* (
					-distancia_radial ** 4
					- 1.5 * planeta.J2 * distancia_radial ** 2 * planeta.raio_equatorial ** 2
					+ 1.875 * planeta.J4 * planeta.raio_equatorial ** 4
					- 6 * planeta.J3 * distancia_radial * planeta.raio_equatorial ** 3 * np.cos(colatitude)
					+ (
								4.5 * planeta.J2 * distancia_radial ** 2 * planeta.raio_equatorial ** 2 - 18.75 * planeta.J4 * planeta.raio_equatorial ** 4) * np.cos(
				colatitude) ** 2
					+ 10 * planeta.J3 * distancia_radial * planeta.raio_equatorial ** 3 * np.cos(
				colatitude) ** 3
					+ 21.875 * planeta.J4 * planeta.raio_equatorial ** 4 * np.cos(colatitude) ** 4
			)
	)

	# Calculando a componente colatitudinal da gravidade
	componente_colatitudinal = (
			(1 / distancia_radial ** 6)
			* 3
			* planeta.gravidade_padrao_nivel_do_mar
			* planeta.constante_gravitacional
			* planeta.raio_equatorial ** 2
			* (
					-0.5 * planeta.J3 * distancia_radial * planeta.raio_equatorial
					+ (
								planeta.J2 * distancia_radial ** 2 - 2.5 * planeta.J4 * planeta.raio_equatorial ** 2) * np.cos(
				colatitude)
					+ 2.5 * planeta.J3 * distancia_radial * planeta.raio_equatorial * np.cos(colatitude) ** 2
					+ 5.83333 * planeta.J4 * planeta.raio_equatorial ** 2 * np.cos(colatitude) ** 3
			)
			* np.sin(colatitude)
	)

	return componente_radial, componente_colatitudinal
