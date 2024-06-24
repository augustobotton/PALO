import numpy as np


def calcular_gravidade_axisimetrico(r, delta, planeta):
	mut = planeta.mut
	Requat = planeta.raio_equatorial
	J2 = planeta.J2
	J3 = planeta.J3
	J4 = planeta.J4

	phi = np.pi / 2 - delta

	gr = -(mut) / (r ** 2) - (1 / r) * mut * (
			-((J2 * Requat ** 2 * (-1 + 3 * np.cos(phi) ** 2)) / (r ** 3)) - (
			1.5 * J3 * Requat ** 3 * (-3 * np.cos(phi) + 5 * np.cos(phi) ** 3)) / (r ** 4) - (
					J4 * Requat ** 4 * (3 - 30 * np.cos(phi) ** 2 + 35 * np.cos(phi) ** 4)) / (
					2 * r ** 5)) + (mut / (r ** 2)) * (
				 (0.5 * J2 * Requat ** 2 * (-1 + 3 * np.cos(phi) ** 2)) / (r ** 2) + (
				 0.5 * J3 * Requat ** 3 * (-3 * np.cos(phi) + 5 * np.cos(phi) ** 3)) / (r ** 3) + (
						 J4 * Requat ** 4 * (3 - 30 * np.cos(phi) ** 2 + 35 * np.cos(phi) ** 4)) / (
						 8 * r ** 4))

	gphi = -(1 / r ** 2) * mut * (-((3 * J2 * Requat ** 2 * np.cos(phi) * np.sin(phi)) / r ** 2) + (
			0.5 * J3 * Requat ** 3 * (3 * np.sin(phi) - 15 * np.cos(phi) ** 2 * np.sin(phi))) / r ** 3 + (
										  J4 * Requat ** 4 * (
										  60 * np.cos(phi) * np.sin(phi) - 140 * np.cos(
									  phi) ** 3 * np.sin(phi))) / (8 * r ** 4))

	gc = -gr
	gdel = -gphi

	return gc, gdel
