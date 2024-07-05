import numpy as np
from src.domain.utilidades_mecanica_orbital.Utilidades.Converte import Vrel2Vine


# Modelo do propulsor de foguete de 1 estágio
def propulsor_1_estagio(t, tempos_de_ignicao, tempos_de_fim_de_queima, tempos_de_separacao, impulso_especifico,
                        massa_propelente_estagios, massa_estrutural_por_estagio
                        , massa_inicial_do_foguete, g):
    # Cálculo da massa e tração em foguete com 1 estágio e carga util
    if t <= tempos_de_ignicao[0]:
        # Antes da ignição
        m = massa_inicial_do_foguete  # Massa inicial
        ft = 0  # Força propulsiva nula
    elif t <= tempos_de_fim_de_queima[0]:
        # Taxa de queima contínua
        md = -massa_propelente_estagios[0] / (tempos_de_fim_de_queima[0] - tempos_de_ignicao[0])
        # Está queimando o primeiro estágio
        m = massa_inicial_do_foguete + md * (t - tempos_de_ignicao[0])
        # Força propulsiva constante
        ft = -g * impulso_especifico[0] * md
    elif t <= tempos_de_separacao[0]:
        # Entre a queima e a separação
        m = massa_inicial_do_foguete - massa_propelente_estagios[0]
        ft = 0
    else:
        # Após a separação do motor foguete
        m = massa_inicial_do_foguete - massa_propelente_estagios[0] - massa_estrutural_por_estagio[0]
        ft = 0
    return ft, m


# Modelo dos propulsores de foguete de 2 estágios
def propulsor_2_estagios(t, tempos_de_ignicao, tempos_de_fim_de_queima, tempos_de_separacao, impulso_especifico,
                         massa_propelente_estagios, massa_estrutural_por_estagio, massa_inicial_do_foguete, g):
    # Cálculo da massa e tração em foguete com 2 estágios e carga util
    if t <= tempos_de_ignicao[0]:
        # Antes da ignição
        m = massa_inicial_do_foguete  # Massa inicial
        ft = 0  # Força propulsiva nula
    elif t <= tempos_de_fim_de_queima[0]:
        # Taxa de queima contínua
        md = -massa_propelente_estagios[0] / (tempos_de_fim_de_queima[0] - tempos_de_ignicao[0])
        # Está queimando o primeiro estágio
        m = massa_inicial_do_foguete + md * (t - tempos_de_ignicao[0])
        # Força propulsiva constante
        ft = -g * impulso_especifico[0] * md
    elif t <= tempos_de_separacao[0]:
        # Entre a queima e a separação
        m = massa_inicial_do_foguete - massa_propelente_estagios[0]
        ft = 0
    elif t <= tempos_de_ignicao[1]:
        # Entre a separação e a ignição
        m = massa_inicial_do_foguete - massa_propelente_estagios[0] - massa_estrutural_por_estagio[0]
        ft = 0
    elif t <= tempos_de_fim_de_queima[1]:
        # Taxa de queima contínua no segundo estágio
        md = -massa_propelente_estagios[1] / (tempos_de_fim_de_queima[1] - tempos_de_ignicao[1])
        # Durante a queima do segundo estágio
        m02 = massa_inicial_do_foguete - massa_propelente_estagios[0] - massa_estrutural_por_estagio[0]
        m = m02 + md * (t - tempos_de_ignicao[1])
        # Força propulsiva constante
        ft = -g * impulso_especifico[1] * md
    elif t <= tempos_de_separacao[1]:
        # Após a queima do segundo estágio e antes da separação do mesmo
        m = massa_inicial_do_foguete - massa_propelente_estagios[0] - massa_estrutural_por_estagio[0] - \
            massa_propelente_estagios[1]
        ft = 0
    else:
        # Após a separação do segundo estágio
        m = massa_inicial_do_foguete - massa_propelente_estagios[0] - massa_estrutural_por_estagio[0] - \
            massa_propelente_estagios[1] - massa_estrutural_por_estagio[1]
        ft = 0
    return ft, m


# Modelo dos propulsores de foguete de 3 estágios
def propulsor_3_estagios(t, tempos_de_ignicao, tempos_de_fim_de_queima, tempos_de_separacao, impulso_especifico,
                         massa_propelente_estagios, massa_estrutural_por_estagio, massa_inicial_do_foguete, g):
    # Cálculo da massa e tração
    if t <= tempos_de_ignicao[0]:
        # Antes da ignição
        m = massa_inicial_do_foguete  # Massa inicial
        ft = 0  # Força propulsiva nula
    elif t <= tempos_de_fim_de_queima[0]:
        # Taxa de queima contínua
        md = -massa_propelente_estagios[0] / (tempos_de_fim_de_queima[0] - tempos_de_ignicao[0])
        # Está queimando o primeiro estágio
        m = massa_inicial_do_foguete + md * (t - tempos_de_ignicao[0])
        # Força propulsiva constante
        ft = -g * impulso_especifico[0] * md
    elif t <= tempos_de_separacao[0]:
        # Entre a queima e a separação
        m = massa_inicial_do_foguete - massa_propelente_estagios[0]
        ft = 0
    elif t <= tempos_de_ignicao[1]:
        # Entre a separação e a ignição
        m = massa_inicial_do_foguete - massa_propelente_estagios[0] - massa_estrutural_por_estagio[0]
        ft = 0
    elif t <= tempos_de_fim_de_queima[1]:
        # Taxa de queima contínua no segundo estágio
        md = -massa_propelente_estagios[1] / (tempos_de_fim_de_queima[1] - tempos_de_ignicao[1])
        # Durante a queima do segundo estágio
        m02 = massa_inicial_do_foguete - massa_propelente_estagios[0] - massa_estrutural_por_estagio[0]
        m = m02 + md * (t - tempos_de_ignicao[1])
        # Força propulsiva constante
        ft = -g * impulso_especifico[1] * md
    elif t <= tempos_de_separacao[1]:
        # Após a queima do segundo estágio e antes da separação do mesmo
        m = massa_inicial_do_foguete - massa_propelente_estagios[0] - massa_estrutural_por_estagio[0] - \
            massa_propelente_estagios[1]
        ft = 0
    elif t <= tempos_de_ignicao[2]:
        # Entre a separação e a ignição
        m = massa_inicial_do_foguete - massa_propelente_estagios[0] - massa_estrutural_por_estagio[0] - \
            massa_propelente_estagios[1] - massa_estrutural_por_estagio[1]
        ft = 0
    elif t <= tempos_de_fim_de_queima[2]:
        # Taxa de queima contínua no terceiro estágio
        md = -massa_propelente_estagios[2] / (tempos_de_fim_de_queima[2] - tempos_de_ignicao[2])
        # Durante a queima do terceiro estágio
        m03 = massa_inicial_do_foguete - massa_propelente_estagios[0] - massa_estrutural_por_estagio[0] - \
              massa_propelente_estagios[1] - massa_estrutural_por_estagio[1]
        m = m03 + md * (t - tempos_de_ignicao[2])
        # Força propulsiva constante
        ft = -g * impulso_especifico[2] * md
    elif t <= tempos_de_separacao[2]:
        # Após a queima do terceiro estágio e antes da separação do mesmo
        m = massa_inicial_do_foguete - massa_propelente_estagios[0] - massa_estrutural_por_estagio[0] - \
            massa_propelente_estagios[1] - massa_estrutural_por_estagio[1] - massa_propelente_estagios[2]
        ft = 0
    else:
        # Após a separação do terceiro estágio
        m = massa_inicial_do_foguete - massa_propelente_estagios[0] - massa_estrutural_por_estagio[0] - \
            massa_propelente_estagios[1] - massa_estrutural_por_estagio[1] - massa_propelente_estagios[2] - \
            massa_estrutural_por_estagio[2]
        ft = 0
    return ft, m


def propulsor_3_estagios_2ig(t, tempos_de_ignicao, tempos_de_fim_de_queima, tempos_de_separacao, impulso_especifico,
                             massa_propelente_estagios, massa_estrutural_por_estagio, massa_inicial_do_foguete,
                             massa_de_carga_util, g):
    # Cálculo da massa e tração
    if t <= tempos_de_ignicao[0]:
        # Antes da ignição
        m = massa_inicial_do_foguete  # Massa inicial
        ft = 0  # Força propulsiva nula
    elif t <= tempos_de_fim_de_queima[0]:
        # Taxa de queima contínua
        md = -massa_propelente_estagios[0] / (tempos_de_fim_de_queima[0] - tempos_de_ignicao[0])
        # Está queimando o primeiro estágio
        m = massa_inicial_do_foguete + md * (t - tempos_de_ignicao[0])
        # Força propulsiva constante
        ft = -g * impulso_especifico[0] * md
    elif t <= tempos_de_separacao[0]:
        # Entre a queima e a separação
        m = massa_inicial_do_foguete - massa_propelente_estagios[0]
        ft = 0
    elif t <= tempos_de_ignicao[1]:
        # Entre a separação e a ignição
        m = massa_inicial_do_foguete - massa_propelente_estagios[0] - massa_estrutural_por_estagio[0]
        ft = 0
    elif t <= tempos_de_fim_de_queima[1]:
        # Taxa de queima contínua no segundo estágio
        md = -massa_propelente_estagios[1] / (tempos_de_fim_de_queima[1] - tempos_de_ignicao[1])
        # Durante a queima do segundo estágio
        m02 = massa_inicial_do_foguete - massa_propelente_estagios[0] - massa_estrutural_por_estagio[0]
        m = m02 + md * (t - tempos_de_ignicao[1])
        # Força propulsiva constante
        ft = -g * impulso_especifico[1] * md
    elif t <= tempos_de_separacao[1]:
        # Após a queima do segundo estágio e antes da separação do mesmo
        m = massa_inicial_do_foguete - massa_propelente_estagios[0] - massa_estrutural_por_estagio[0] - \
            massa_propelente_estagios[1]
        ft = 0
    elif t <= tempos_de_ignicao[2]:
        # Entre a separação e a ignição
        m = massa_inicial_do_foguete - massa_propelente_estagios[0] - massa_estrutural_por_estagio[0] - \
            massa_propelente_estagios[1] - massa_estrutural_por_estagio[1]
        ft = 0
    elif t <= tempos_de_fim_de_queima[2]:
        # Taxa de queima contínua no terceiro estágio - primeira ignição
        md = -massa_propelente_estagios[2] / (tempos_de_fim_de_queima[2] - tempos_de_ignicao[2])
        # Durante a queima do terceiro estágio
        m03 = massa_inicial_do_foguete - massa_propelente_estagios[0] - massa_estrutural_por_estagio[0] - \
              massa_propelente_estagios[1] - massa_estrutural_por_estagio[1]
        m = m03 + md * (t - tempos_de_ignicao[2])
        # Força propulsiva constante
        ft = -g * impulso_especifico[2] * md
    elif t <= tempos_de_ignicao[3]:
        # Antes da nova queima do terceiro estágio
        m = massa_inicial_do_foguete - massa_propelente_estagios[0] - massa_estrutural_por_estagio[0] - \
            massa_propelente_estagios[1] - massa_estrutural_por_estagio[1] - massa_propelente_estagios[2]
        ft = 0
    elif t <= tempos_de_fim_de_queima[3]:
        # Taxa de queima contínua no terceiro estágio - segunda ignição
        md = -massa_propelente_estagios[3] / (tempos_de_fim_de_queima[3] - tempos_de_ignicao[3])
        # Durante a queima do terceiro estágio
        m03 = massa_inicial_do_foguete - massa_propelente_estagios[0] - massa_estrutural_por_estagio[0] - \
              massa_propelente_estagios[1] - massa_estrutural_por_estagio[1] - massa_propelente_estagios[2]
        m = m03 + md * (t - tempos_de_ignicao[3])
        # Força propulsiva constante
        ft = -g * impulso_especifico[2] * md
    elif t <= tempos_de_separacao[2]:
        # Após a queima do terceiro estágio e antes da separação do mesmo
        m = massa_inicial_do_foguete - massa_propelente_estagios[0] - massa_estrutural_por_estagio[0] - \
            massa_propelente_estagios[1] - massa_estrutural_por_estagio[1] - massa_propelente_estagios[2] - \
            massa_propelente_estagios[3]
        ft = 0
    else:
        # Após a separação do terceiro estágio
        m = massa_de_carga_util
        ft = 0
    return ft, m
