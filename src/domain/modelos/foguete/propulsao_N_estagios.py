import numpy as np
import parametros

from src.domain.orbitalUtils.Vrel2Vine import Vrel2Vine


def propulsao_N_estagios(tempo, vetor_de_estados):
    ti = parametros.ti
    tq = parametros.tq
    ts = parametros.tempo_limite_separacao
    Isp = parametros.impulso_especifico_por_estagio
    mp = parametros.mp
    ms = parametros.massa_estrutural_por_estagio
    g = parametros.gravidade_padrao_nivel_do_mar
    m0 = parametros.m0
    Re = parametros.raio_equatorial
    we = parametros.velocidade_inercial_de_rotação_da_terra
    mL = parametros.massa_de_carga_util
    # Função para cálculo dos parâmetros propulsivos em função do tempo
    # Veículo de até 3 estágios
    # Entrada
    # t (s): tempo
    # Saídas
    # ft (N): força propulsiva
    # m (kg): massa do foguete em função do tempo
    # Hipoteses
    # - Em cada estagio, eh assumida taxa de queima continua, ou seja, nao ha
    # controle da queima e a mesma eh assumida uniforme do inicio ao fim do
    # propelente
    # - A tracao de cada estagio eh assumida como um pulso retangular perfeito,
    # ou seja, quando acionado, o propulsor vai de tracao zero ateh a maxima,
    # permanecendo nesse patamar constante. Ao fim da queima, a tracao cai
    # instanteamente a zero

    # As variáveis propulsivas dependem do estagio atual
    # Numero de estagios
    N = len(ti)

    if N == 1:
        ft, m = propulsor_1_estagio(tempo, ti, tq, ts, Isp, mp, ms, m0, g)
    elif N == 2:
        ft, m = propulsor_2_estagios(tempo, ti, tq, ts, Isp, mp, ms, m0, g)
    elif N == 3:
        ft, m = propulsor_3_estagios(tempo, ti, tq, ts, Isp, mp, ms, m0, g)
    elif N == 4:
        ft, m = propulsor_3_estagios_2ig(tempo, ti, tq, ts, Isp, mp, ms, m0, mL, g)

    V = vetor_de_estados[0]
    A = vetor_de_estados[1]
    phi = vetor_de_estados[2]
    r = vetor_de_estados[3]
    delta = vetor_de_estados[4]

    # Altitude
    h = r - Re
    if h < 200e3:
        epsl = 0
        mu = 0
    else:
        _, phii, Ai = Vrel2Vine(V, phi, A, we, r, delta)

        mu = np.arcsin(np.cos(A) * np.cos(phii) * np.sin(Ai) - np.sin(A) * np.cos(phii) * np.cos(Ai))
        epsl = -np.arctan2(-np.cos(phi) * np.sin(phii) + np.sin(phi) * np.sin(A) * np.cos(phii) * np.sin(Ai) +
                           np.sin(phi) * np.cos(A) * np.cos(phii) * np.cos(Ai), np.sin(phi) * np.sin(phii) +
                           np.cos(phi) * np.sin(A) * np.cos(phii) * np.sin(Ai) + np.cos(phi) * np.cos(A) * np.cos(
            phii) * np.cos(Ai))

    return ft, m, mu, epsl


def propulsor_1_estagio(t, ti, tq, ts, Isp, mp, ms, m0, g):
    # Modelo do propulsor de foguete de 1 estágio
    # Cálculo da massa e tração em foguete com 1 estagio e carga util

    if t <= ti[0]:
        # Antes da ignição
        m = m0  # Massa inicial
        ft = 0  # Força propulsiva nula
    elif t <= tq[0]:
        # Taxa de queima contínua
        md = -mp[0] / (tq[0] - ti[0])
        # Está queimando o primeiro estágio
        m = m0 + md * (t - ti[0])
        # Força propulsiva constante
        ft = -g * Isp[0] * md
    elif t <= ts[0]:
        # Entre a queima e a separação
        m = m0 - mp[0]
        ft = 0
    else:
        # Após a separação do motor foguete
        m = m0 - mp[0] - ms[0]
        ft = 0

    return ft, m


def propulsor_2_estagios(t, ti, tq, ts, Isp, mp, ms, m0, g):
    # Modelo dos propulsores de foguete de 2 estágios
    # Cálculo da massa e tração em foguete com 2 estagios e carga util

    if t <= ti[0]:
        # Antes da ignição
        m = m0  # Massa inicial
        ft = 0  # Força propulsiva nula
    elif t <= tq[0]:
        # Taxa de queima contínua
        md = -mp[0] / (tq[0] - ti[0])
        # Está queimando o primeiro estágio
        m = m0 + md * (t - ti[0])
        # Força propulsiva constante
        ft = -g * Isp[0] * md
    elif t <= ts[0]:
        # Entre a queima e a separação
        m = m0 - mp[0]
        ft = 0
    elif t <= ti[1]:
        # Entre a separação e a ignição
        m = m0 - mp[0] - ms[0]
        ft = 0
    elif t <= tq[1]:
        # Taxa de queima contínua no segundo estágio
        md = -mp[1] / (tq[1] - ti[1])
        # Durante a queima do segundo estágio
        m02 = m0 - mp[0] - ms[0]
        m = m02 + md * (t - ti[1])
        # Força propulsiva constante
        ft = -g * Isp[1] * md
    elif t <= ts[1]:
        # Após a queima do segundo estágio e antes da separação do mesmo
        m = m0 - mp[0] - ms[0] - mp[1]
        ft = 0
    else:
        # Após a separação do segundo estágio
        m = m0 - mp[0] - ms[0] - mp[1] - ms[1]
        ft = 0

    return ft, m


def propulsor_3_estagios(t, ti, tq, ts, Isp, mp, ms, m0, g):
    # Modelo dos propulsores de foguete de 3 estágios
    # Cálculo da massa e tração

    if t <= ti[0]:
        # Antes da ignição
        m = m0  # Massa inicial
        forca_propulsiva = 0  # Força propulsiva nula
    elif t <= tq[0]:
        # Taxa de queima contínua
        md = -mp[0] / (tq[0] - ti[0])
        # Está queimando o primeiro estágio
        m = m0 + md * (t - ti[0])
        # Força propulsiva constante
        forca_propulsiva = -g * Isp[0] * md
    elif t <= ts[0]:
        # Entre a queima e a separação
        m = m0 - mp[0]
        forca_propulsiva = 0
    elif t <= ti[1]:
        # Entre a separação e a ignição
        m = m0 - mp[0] - ms[0]
        forca_propulsiva = 0
    elif t <= tq[1]:
        # Taxa de queima contínua no segundo estágio
        md = -mp[1] / (tq[1] - ti[1])
        # Durante a queima do segundo estágio
        m02 = m0 - mp[0] - ms[0]
        m = m02 + md * (t - ti[1])
        # Força propulsiva constante
        forca_propulsiva = -g * Isp[1] * md
    elif t <= ts[1]:
        # Após a queima do segundo estágio e antes da separação do mesmo
        m = m0 - mp[0] - ms[0] - mp[1]
        forca_propulsiva = 0
    elif t <= ti[2]:
        # Entre a separação e a ignição
        m = m0 - mp[0] - ms[0] - mp[1] - ms[1]
        forca_propulsiva = 0
    elif t <= tq[2]:
        # Taxa de queima contínua no terceiro estágio
        md = -mp[2] / (tq[2] - ti[2])
        # Durante a queima do terceiro estágio
        m03 = m0 - mp[0] - ms[0] - mp[1] - ms[1]
        m = m03 + md * (t - ti[2])
        # Força propulsiva constante
        forca_propulsiva = -g * Isp[2] * md
    elif t <= ts[2]:
        # Após a queima do terceiro estágio e antes da separação do mesmo
        m = m0 - mp[0] - ms[0] - mp[1] - ms[1] - mp[2]
        forca_propulsiva = 0
    else:
        # Após a separação do terceiro estágio
        m = m0 - mp[0] - ms[0] - mp[1] - ms[1] - mp[2] - ms[2]
        forca_propulsiva = 0

    return forca_propulsiva, m


def propulsor_3_estagios_2ig(t, ti, tq, ts, Isp, mp, ms, m0, mL, g):
    if t <= ti[0]:
        # Antes da ignição
        massa_em_funcao_do_tempo = m0  # Massa inicial
        forca_propulsiva = 0  # Força propulsiva nula
    elif t <= tq[0]:
        # Taxa de queima contínua
        md = -mp[0] / (tq[0] - ti[0])
        # Está queimando o primeiro estágio
        massa_em_funcao_do_tempo = m0 + md * (t - ti[0])
        # Força propulsiva constante
        forca_propulsiva = -g * Isp[0] * md
    elif t <= ts[0]:
        # Entre a queima e a separação
        massa_em_funcao_do_tempo = m0 - mp[0]
        forca_propulsiva = 0
    elif t <= ti[1]:
        # Entre a separação e a ignição
        massa_em_funcao_do_tempo = m0 - mp[0] - ms[0]
        forca_propulsiva = 0
    elif t <= tq[1]:
        # Taxa de queima contínua no segundo estágio
        md = -mp[1] / (tq[1] - ti[1])
        # Durante a queima do segundo estágio
        m02 = m0 - mp[0] - ms[0]
        massa_em_funcao_do_tempo = m02 + md * (t - ti[1])
        # Força propulsiva constante
        forca_propulsiva = -g * Isp[1] * md
    elif t <= ts[1]:
        # Após a queima do segundo estágio e antes da separação do mesmo
        massa_em_funcao_do_tempo = m0 - mp[0] - ms[0] - mp[1]
        forca_propulsiva = 0
    elif t <= ti[2]:
        # Entre a separação e a ignição
        massa_em_funcao_do_tempo = m0 - mp[0] - ms[0] - mp[1] - ms[1]
        forca_propulsiva = 0
    elif t <= tq[2]:
        # Taxa de queima contínua no terceiro estágio - primeira ignicao
        md = -mp[2] / (tq[2] - ti[2])
        # Durante a queima do terceiro estágio
        m03 = m0 - mp[0] - ms[0] - mp[1] - ms[1]
        massa_em_funcao_do_tempo = m03 + md * (t - ti[2])
        # Força propulsiva constante
        forca_propulsiva = -g * Isp[2] * md
    elif t <= ti[3]:
        # Antes da nova queima do terceiro estagio
        massa_em_funcao_do_tempo = m0 - mp[0] - ms[0] - mp[1] - ms[1] - mp[2]
        forca_propulsiva = 0
    elif t <= tq[3]:
        # Taxa de queima contínua no terceiro estágio - segunda ignicao
        md = -mp[3] / (tq[3] - ti[3])
        # Durante a queima do terceiro estágio
        m03 = m0 - mp[0] - ms[0] - mp[1] - ms[1] - mp[2]
        massa_em_funcao_do_tempo = m03 + md * (t - ti[3])
        # Força propulsiva constante
        forca_propulsiva = -g * Isp[2] * md
    elif t <= ts[2]:
        # Após a queima do terceiro estágio e antes da separação do mesmo
        massa_em_funcao_do_tempo = m0 - mp[0] - ms[0] - mp[1] - ms[1] - mp[2] - mp[3]
        forca_propulsiva = 0
    else:
        # Após a separação do terceiro estágio
        massa_em_funcao_do_tempo = mL
        forca_propulsiva = 0

    return forca_propulsiva, massa_em_funcao_do_tempo
