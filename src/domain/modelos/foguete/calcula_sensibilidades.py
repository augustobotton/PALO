import numpy as np


def calcula_sensibilidades(massa_inicial, massa_propelente, massa_carga_util, velocidade_exaustao):
    """
    Calcula as sensibilidades das massas estruturais e de propelente em relação à carga útil.

    Args:
        massa_inicial (np.array): Massas no início da queima de cada estágio.
        massa_propelente (np.array): Massas de propelente consumidas em cada estágio.
        massa_carga_util (float): Massa da carga útil.
        velocidade_exaustao (np.array): Velocidade de exaustão de cada estágio.

    Resultados:
        Imprime as massas estruturais, de propelente, finais e as derivadas calculadas.
    """
    # Determinação das massas estruturais de cada estágio
    num_estagios = len(massa_inicial)
    massa_estrutural = np.zeros(num_estagios)
    massa_estrutural[-1] = massa_inicial[-1] - massa_carga_util - massa_propelente[-1]
    massa_estrutural[:-1] = massa_inicial[:-1] - massa_inicial[1:] - massa_propelente[:-1]

    # Massas ao final da queima de cada estágio
    massa_final = massa_inicial - massa_propelente

    # Mostrar as massas
    print(f"Massa da carga útil (kg): {massa_carga_util}")
    print(f"Massas no início da queima de cada estágio: {massa_inicial}")
    print(f"Massas no final da queima de cada estágio: {massa_final}")
    print(f"Massas estruturais de cada estágio: {massa_estrutural}")
    print(f"Massas de propelente de cada estágio: {massa_propelente}")

    # Cálculo das derivadas de massa de carga útil com respeito às massas estruturais
    derivadas_massa_estrutural = calcular_derivadas_massa_estrutural(massa_propelente, massa_estrutural,
                                                                     massa_carga_util, velocidade_exaustao)
    print(f"Vetor de derivadas da carga útil com respeito às massas estruturais: {derivadas_massa_estrutural}")

    # Cálculo das derivadas de massa de carga útil com respeito às massas de propelente
    derivadas_massa_propelente = calcular_derivadas_massa_propelente(massa_propelente, massa_estrutural,
                                                                     massa_carga_util, velocidade_exaustao)
    print(f"Vetor de derivadas da carga útil com respeito às massas de propelente: {derivadas_massa_propelente}")


def calcular_derivadas_massa_estrutural(massa_propelente, massa_estrutural, massa_carga_util, velocidade_exaustao):
    """
    Calcula as derivadas da massa de carga útil com respeito às massas estruturais.

    Args:
        massa_propelente (np.array): Massas de propelente consumidas em cada estágio.
        massa_estrutural (np.array): Massas estruturais de cada estágio.
        massa_carga_util (float): Massa da carga útil.
        velocidade_exaustao (np.array): Velocidade de exaustão de cada estágio.

    Returns:
        np.array: Derivadas da massa de carga útil com respeito às massas estruturais.
    """
    num_estagios = len(massa_propelente)
    massa_final = np.zeros(num_estagios)
    massa_inicial = np.zeros(num_estagios)
    for i in range(num_estagios):
        massa_inicial[i] = massa_carga_util + np.sum(massa_propelente[i:]) + np.sum(massa_estrutural[i:])
        massa_final[i] = massa_inicial[i] - massa_propelente[i]

    derivadas = np.zeros(num_estagios)
    for k in range(num_estagios):
        derivadas[k] = -np.sum(
            velocidade_exaustao[:k + 1] * (1.0 / massa_inicial[:k + 1] - 1.0 / massa_final[:k + 1])) / np.sum(
            velocidade_exaustao * (1.0 / massa_inicial - 1.0 / massa_final))

    return derivadas


def calcular_derivadas_massa_propelente(massa_propelente, massa_estrutural, massa_carga_util, velocidade_exaustao):
    """
    Calcula as derivadas da massa de carga útil com respeito às massas de propelente.

    Args:
        massa_propelente (np.array): Massas de propelente consumidas em cada estágio.
        massa_estrutural (np.array): Massas estruturais de cada estágio.
        massa_carga_util (float): Massa da carga útil.
        velocidade_exaustao (np.array): Velocidade de exaustão de cada estágio.

    Returns:
        np.array: Derivadas da massa de carga útil com respeito às massas de propelente.
    """
    num_estagios = len(massa_propelente)
    massa_final = np.zeros(num_estagios)
    massa_inicial = np.zeros(num_estagios)
    for i in range(num_estagios):
        massa_inicial[i] = massa_carga_util + np.sum(massa_propelente[i:]) + np.sum(massa_estrutural[i:])
        massa_final[i] = massa_inicial[i] - massa_propelente[i]

    derivadas = np.zeros(num_estagios)
    for k in range(num_estagios):
        if k == 0:
            derivadas[k] = -velocidade_exaustao[k] / massa_inicial[k] / np.sum(
                velocidade_exaustao * (1.0 / massa_inicial - 1.0 / massa_final))
        else:
            derivadas[k] = -(velocidade_exaustao[k] / massa_inicial[k] + np.sum(
                velocidade_exaustao[:k] * (1.0 / massa_inicial[:k] - 1.0 / massa_final[:k]))) / np.sum(
                velocidade_exaustao * (1.0 / massa_inicial - 1.0 / massa_final))

    return derivadas


# Dados de entrada do exemplo 8.4
massa_inicial = np.array([109426.6769684211, 52847.72960000001, 19806.1060, 6046.9099])
massa_propelente = np.array([55000, 28978.709948, 12796.05244274179, 4568.81457])
massa_carga_util = 1124.811598501497
velocidade_exaustao = 1000 * np.array([2.363318181818182, 2.8449, 2.8449, 4.463550000000001])

# Executa a função principal para calcular as sensibilidades com os dados de entrada fornecidos
calcula_sensibilidades(massa_inicial, massa_propelente, massa_carga_util, velocidade_exaustao)
