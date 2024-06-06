import numpy as np

def det_orbita(tempo_observacao, posicao_relativa_ao_primario_referencial_celeste, velocidade_relativa_ao_primario_referencial_celeste, parametro_gravitacional_corpo_primario):
    # Função para determinar parâmetros orbitais semi_eixo_maior partir de uma observação de
    # posição excentricidade_orbita outra de velocidade, sendo as mesmas tomadas em relação ao
    # primário em um problema de dois corpos excentricidade_orbita escritas no referencial
    # celestial
    # Entradas
    # t0: tempo em que semi_eixo_maior observacao foi feita (s)
    # rc0: vetor posição em relação ao primário escrito no referencial
    # celeste (m ou km)
    # vc0: vetor velocidade em relação ao primário escrito no referencial
    # celeste (m/s ou km/s - unidades coerentes com o vetor posição). Deve
    # ser tomado no mesmo instante da medida da posição.
    # mu: parâmetro gravitacional padrão do corpo primário (m^3/s^2 ou km^3/s^2
    # - unidades coerentes com semi_eixo_maior posição excentricidade_orbita velocidade)
    # Saídas
    # par_orb: vetor de parâmetros orbitais
    # semi_eixo_maior=par_orb[0]: semi eixo maior da órbita (m ou km - depende das unidades
    # de entrada)
    # excentricidade_orbita=par_orb[1]: excentricidade da órbita (adimensional)
    # tempo_de_perigeu=par_orb[2]: tempo de perigeu da órbita (segundos)
    # ascenção_reta_do_nodo_ascendente=par_orb[3]: ascenção reta do nodo ascendente (rad). Direção da
    # linha dos nodos da órbita dada em relação ao eixo X no plano XY do
    # referencial celeste
    # i=par_orb[4]: inclinação (rad). Inclinação da órbita dada em relação ao
    # plano XY do referencial celeste.
    # argumento_de_perigeu=par_orb[5]: argumento de perigeu (rad). Em relação à linha dos
    # nodos, medido no plano orbital.

    ## Cálculos
    # Distância radial ao primário no instante observado
    posicao_relativa_ao_primario_referencial_celeste = posicao_relativa_ao_primario_referencial_celeste.T  # para funcionar no np.cross
    velocidade_relativa_ao_primario_referencial_celeste = velocidade_relativa_ao_primario_referencial_celeste.T
    posicao_relativa_ao_primario_normalizada = np.linalg.norm(posicao_relativa_ao_primario_referencial_celeste)

    quantidade_de_movimento_espceifica_ref_celeste = np.cross(posicao_relativa_ao_primario_referencial_celeste, velocidade_relativa_ao_primario_referencial_celeste)  # como rc0 excentricidade_orbita vc0 foram transpostos, quantidade_de_movimento_espceifica_ref_celeste vai ser 1x3 tambem

    excentricidade_sistema_celeste = np.cross(velocidade_relativa_ao_primario_referencial_celeste, quantidade_de_movimento_espceifica_ref_celeste) / parametro_gravitacional_corpo_primario - posicao_relativa_ao_primario_referencial_celeste / posicao_relativa_ao_primario_normalizada  # voltando pra (3x1)
    excentricidade_orbita = np.linalg.norm(excentricidade_sistema_celeste)
    # Módulo do vetor quantidade_de_movimento_espceifica_ref_celeste
    modulo_quatidade_de_movimento_especifica_no_referencial_celeste = np.linalg.norm(quantidade_de_movimento_espceifica_ref_celeste)

    # Parâmetro da órbita dada
    parametro_orbita_dada = modulo_quatidade_de_movimento_especifica_no_referencial_celeste ** 2 / parametro_gravitacional_corpo_primario

    # Semi eixo maior
    semi_eixo_maior = parametro_orbita_dada / (1 - excentricidade_orbita ** 2)

    # Vetor parâmetro no referencial celeste
    parametro_referencial_celeste = parametro_orbita_dada * np.cross(quantidade_de_movimento_espceifica_ref_celeste, excentricidade_sistema_celeste) / (modulo_quatidade_de_movimento_especifica_no_referencial_celeste * excentricidade_orbita)

    # Anomalia verdadeira
    costheta = (parametro_orbita_dada - posicao_relativa_ao_primario_normalizada) / (excentricidade_orbita * posicao_relativa_ao_primario_normalizada)
    sintheta = np.dot(np.squeeze(posicao_relativa_ao_primario_referencial_celeste), np.squeeze(parametro_referencial_celeste)) / (posicao_relativa_ao_primario_normalizada * parametro_orbita_dada)
    anomalia_verdadeira = np.arctan2(sintheta, costheta)

    # O tempo de perigeu depende do tipo de órbita
    if (0 <= excentricidade_orbita) and (excentricidade_orbita < 1):
        tipo = 'excentricidade_orbita' # Órbita elíptica
    elif excentricidade_orbita == 1:
        tipo = 'parametro_orbita_dada'  # Órbita parabólica
    else:
        tipo = 'h'  # Órbita hiperbólica

    # Tempo de perigeu
    if tipo == 'excentricidade_orbita':  # Órbita elíptica
        # Movimento médio
        movimento_medio_hiperbolico = np.sqrt(parametro_gravitacional_corpo_primario / semi_eixo_maior ** 3)
        # Anomalia excêntrica
        anomalia_excentrica = 2 * np.arctan(np.sqrt((1 - excentricidade_orbita) / (1 + excentricidade_orbita)) * np.tan(anomalia_verdadeira / 2))
        tempo_de_perigeu = tempo_observacao - (anomalia_excentrica - excentricidade_orbita * np.sin(anomalia_excentrica)) / movimento_medio_hiperbolico
    elif tipo == 'parametro_orbita_dada':  # Órbita parabólica

        tempo_de_perigeu = -((np.tan(anomalia_verdadeira / 2)) ** 3 + 3 * np.tan(anomalia_verdadeira / 2)) / (parametro_gravitacional_corpo_primario / parametro_orbita_dada ** 3) ** (1 / 6)
    else:  # Órbita hiperbólica

        # Movimento médio hiperbólico
        movimento_medio_hiperbolico = np.sqrt(-parametro_gravitacional_corpo_primario / semi_eixo_maior ** 3)
        # Anomalia hiperbólica
        anomalia_hiperbolica = 2 * np.arctanh(np.sqrt((excentricidade_orbita - 1) / (1 + excentricidade_orbita)) * np.tan(
            anomalia_verdadeira / 2))
        tempo_de_perigeu = -(excentricidade_orbita * np.sinh(anomalia_hiperbolica) - anomalia_hiperbolica) / movimento_medio_hiperbolico

    # Linha dos nodos
    # Vetor unitário ao longo do vetor h (no sistema celeste)
    ih = quantidade_de_movimento_espceifica_ref_celeste / modulo_quatidade_de_movimento_especifica_no_referencial_celeste;
    # Vetor unitário ao longo da linha dos nodos (no sistema celeste)
    Kc = np.array([[0], [0], [1]])
    nc = np.cross(Kc.T, ih) / np.linalg.norm(np.cross(Kc.T, ih))
    # Ascenção reta do nodo ascendente
    ascenção_reta_do_nodo_ascendente = np.arctan2(nc[0][1], nc[0][0])
    # Inclinação
    inclinacao = np.arccos(float(np.dot(ih, Kc)))

    # Vetor unitário ao longo do vetor excentricidade (no referencial
    # celeste)
    ie = excentricidade_sistema_celeste / excentricidade_orbita
    # Argumento de perigeu
    cosomega = np.dot(np.squeeze(ie), np.squeeze(nc))
    sinomega = np.dot(np.squeeze(ih), np.squeeze(np.cross(nc, ie)))
    argumento_de_perigeu = np.arctan2(sinomega, cosomega)

    ## Vetor de parâmetros de saída
    parametros_orbitais = [semi_eixo_maior, excentricidade_orbita, tempo_de_perigeu, ascenção_reta_do_nodo_ascendente, inclinacao, argumento_de_perigeu]

    return parametros_orbitais
