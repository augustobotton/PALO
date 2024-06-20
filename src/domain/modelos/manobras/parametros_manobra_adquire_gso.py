import numpy as np

from src.domain.utilidades_mecanica_orbital.orbitalUtils.Converte import Vrel2Vine


def parametros_manobra_adquire_gso(t, m, X, orbita_alvo, foguete):
    """
    Função para calcular os parâmetros da manobra mono impulsiva de aquisição de órbita GSO.
    Deve ser chamada no final da função de dinâmica, pois precisa rodar ao longo do tempo da simulação,
    verificando quando ocorre o apogeu da órbita GTO e determinando os parâmetros para imprimir o impulso
    de velocidade de circularização. Esta função não retorna saídas, mas atualiza parâmetros por variáveis globais,
    que serão usados na função de propulsão.

    Parâmetros:
    - t (s): tempo
    - m (kg): massa do foguete
    - X: vetor de estado [V, A, phi, r, delta]

    """
    we = velocidade_inercial_de_rotação_da_terra
    g = parametros.gravidade_padrao_nivel_do_mar
    velocidade_orbital_desejada = orbita_alvo.velocidade

    sinalPhii = np.sign
    achouApogeu = 0

    ti = parametros.ti
    tq = parametros.tq
    ts = parametros.tempo_limite_separacao
    Tq3 = parametros.Tq3
    Tq31 = parametros.tempo_da_primeira_queima_3_estagio
    Ts3 = parametros.tempo_de_espera_separacao_3_4
    Isp = parametros.impulso_especifico_por_estagio
    mp3 = parametros.mp3

    # Desmenbra o vetor de estado
    V = X[0]
    A = X[1]
    phi = X[2]
    r = X[3]
    delta = X[4]

    # Vetor velocidade inercial
    Vi, phii, _ = Vrel2Vine(V, phi, A, we, r, delta)
    apogeu_orbita_desejada = orbita_alvo.apogeu
    # Realização de uma sequência de testes para verificar a ocorrência do apogeu da órbita GTO.
    # Quando ele ocorre, determina os parâmetros da manobra.
    if r > 0.9 * agso:
        if not achouApogeu:
            if np.sign(phii) != sinalPhii:
                # Se o sinal for diferente, phii passou por zero, o foguete chegou no apogeu
                achouApogeu = 1  # Encontrou o apogeu, ao setar essa variável, só vai entrar aqui uma vez
                parametros.achouApogeu = achouApogeu
                ti[3] = t  # Tempo da segunda ignição do motor do terceiro estágio
                parametros.ti[3] = ti[3]
                parametros.ti = ti
                # Calcula o tempo de duração da queima para propiciar o DeltaV necessário

                DVgso = vgso - Vi  # Quando se nos passa nos testes acima, "Vi" é a velocidade de apogeu da GTO
                parametros.DVgso = DVgso
                mp32 = (m * np.exp(DVgso / (Isp[2] * g)) - m) / np.exp(
                    DVgso / (Isp[2] * g))  # Massa de propelente necessária
                parametros.mp32 = mp32
                Tq32 = (Tq3 * mp32) / mp3  # Duração da queima necessária
                parametros.tempo_da_segunda_queima_3_estagio = Tq32
                # Verifica se o tempo é maior que o máximo, se ocorrer, corrige
                if Tq31 + Tq32 > Tq3:
                    Tq32 = Tq3 - Tq31
                    parametros.tempo_da_segunda_queima_3_estagio = Tq32
                tq[3] = ti[3] + Tq32  # Tempo de fim de queima
                parametros.tq[3] = tq[3]
                parametros.tq = tq
                ts[2] = tq[3] + Ts3  # Tempo de separação
                parametros.tempo_limite_separacao[2] = ts[2]
                parametros.tempo_limite_separacao = ts

    sinalPhii = np.sign(phii)  # Guarda o sinal de phi inercial para verificar mudança na próxima iteração
    parametros.sinalPhii = sinalPhii
