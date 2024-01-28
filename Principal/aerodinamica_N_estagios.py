import modelo_aerodinamico as ma
import parametros


def aerodinamica_multiplos_estagios(tempo, velocidade, altitude, numero_de_mach, numero_knudsen, temperatura,
                                    densidade_do_ar,
                                    constante_do_gas_ideal):
    """
    Calcula as forças aerodinâmicas para foguetes de múltiplos estágios.

    """

    tempo_limite_separacao = parametros.tempo_limite_separacao
    area_de_referencia = parametros.area_de_referencia
    fator_correcao_arrasto = parametros.fator_correcao_arrasto

    # Calcular o coeficiente de arrasto
    coeficiente_arrasto = ma.ModeloAerodinamico(velocidade, altitude, numero_de_mach, numero_knudsen, temperatura,
                                                constante_do_gas_ideal).calcula()
    coeficiente_arrasto_ajustado = fator_correcao_arrasto * coeficiente_arrasto

    # Determinar a área de referência com base no estágio
    area_do_estagio = calcula_area_referencia(tempo, tempo_limite_separacao, area_de_referencia)

    # Calcular as forças
    arrasto = 0.5 * densidade_do_ar * velocidade ** 2 * area_do_estagio * coeficiente_arrasto_ajustado
    forca_sustentacao_lateral = 0
    forca_sustentacao = 0

    return arrasto, forca_sustentacao_lateral, forca_sustentacao


def calcula_area_referencia(tempo, tempo_limite_separacao, area_de_referencia):
    """
    Calcula a área de referência com base no estágio do foguete.
    """
    N = len(tempo_limite_separacao)
    for i in range(N):
        if tempo <= tempo_limite_separacao[i]:
            return area_de_referencia[i]
    return area_de_referencia[-1]  # Padrão para a área do último estágio
