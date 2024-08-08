import numpy as np

"""
Universidade: Universidade Federal de Santa Maria
Curso: Engenharia Aeroespacial
Projeto:  DESENVOLVIMENTO DE UMABIBLIOTECA PYTHON PARA CÁLCULOS DE MECÂNICA ORBITAL E SIMULAÇÃO DE VOO ASCENDENTE
Autor: Augusto Botton Pozzebon
Orientador: Prof. André Luis da Silva
Data: 2024/1
Baseado em: "Orbital Mechanics for Engineering Students" de Howard D. Curtis 3ed
Equações 3.52 e 3.53 
Informações de Contato:
GitHub: https://github.com/augustobotton/PALO
"""



def stumpC(z):
    """
    Calcula a função Stumpff C(z).

    Parâmetros:
        z (float): Valor de entrada para a função Stumpff C.

    Retorna:
        float: Resultado da função Stumpff C(z).
    """

    if z > 0:
        c = (1 - np.cos(np.sqrt(z))) / z
    elif z < 0:
        c = (np.cosh(np.sqrt(-z)) - 1) / (-z)
    else:
        c = 1 / 2

    return c


def stumpS(z):
    """
    Calcula a função Stumpff S(z).

    Parâmetros:
        z (float): Valor de entrada para a função Stumpff S.

    Retorna:
        float: Resultado da função Stumpff S(z).
    """
    if z > 0:
        return (np.sqrt(z) - np.sin(np.sqrt(z))) / (np.sqrt(z)**3)
    elif z < 0:
        return (np.sinh(np.sqrt(-z)) - np.sqrt(-z)) / (np.sqrt(-z)**3)
    else:
        return 1 / 6
