import numpy as np

"""
Universidade: Universidade Federal de Santa Maria
Curso: Engenharia Aeroespacial
Projeto:  DESENVOLVIMENTO DE UMABIBLIOTECA PYTHON PARA CÁLCULOS DE MECÂNICA ORBITAL E SIMULAÇÃO DE VOO ASCENDENTE
Autor: Augusto Botton Pozzebon
Orientador: Prof. André Luis da Silva
Data: 2024/1
Baseado em: "Orbital Mechanics for Engineering Students" de Howard D. Curtis 3ed
Número do Algoritmo: 4.5 (Conversão de Elementos Orbitais Clássicos para Vetor de Estado)

Informações de Contato:

GitHub: https://github.com/augustobotton/PALO
"""
def coe_para_vetor_de_estados(coe, mu):
    """
    This function computes the state vector (r, v) from the classical orbital elements (coe).
    """

    e = coe[1]
    incl = coe[2]
    RA = coe[3]
    w = coe[4]
    TA = coe[5]
    h = coe[6]

    # Position and velocity in the perifocal coordinate system
    rp = (h ** 2 / mu) * (1 / (1 + e * np.cos(TA))) * (
            np.cos(TA) * np.array([1, 0, 0]) + np.sin(TA) * np.array([0, 1, 0]))
    vp = (mu / h) * (-np.sin(TA) * np.array([1, 0, 0]) + (e + np.cos(TA)) * np.array([0, 1, 0]))

    # Rotation matrix about the z-axis through the angle RA
    R3_RA = np.array([
        [np.cos(RA), np.sin(RA), 0],
        [-np.sin(RA), np.cos(RA), 0],
        [0, 0, 1]
    ])

    # Rotation matrix about the x-axis through the angle incl
    R1_i = np.array([
        [1, 0, 0],
        [0, np.cos(incl), np.sin(incl)],
        [0, -np.sin(incl), np.cos(incl)]
    ])

    # Rotation matrix about the z-axis through the angle w
    R3_w = np.array([
        [np.cos(w), np.sin(w), 0],
        [-np.sin(w), np.cos(w), 0],
        [0, 0, 1]
    ])

    # Transformation matrix from perifocal to geocentric equatorial frame
    Q_pX = (R3_w @ R1_i @ R3_RA).T

    # Position and velocity in the geocentric equatorial frame
    r = Q_pX @ rp
    v = Q_pX @ vp

    return r, v