import numpy as np
from math import asin, atan2
from Vrel2Vine import Vrel2Vine
import aquisicao_orbita_gso_na_mao as ent

def propulsao_N_estagios(t, V, A, phi, r, delta):

    def propulsor_1_estagio(t, ti, tq, ts, Isp, mp, ms, m0, g):
        if t <= ti[0]:
            m = m0  # Massa inicial
            ft = 0  # Forca propulsiva nula
        elif t <= tq[0]:
            md = -mp[0] / (tq[0] - ti[0])  # Taxa de queima continua
            m = m0 + md * (t - ti[0])  # Esta queimando o primeiro estagio
            ft = -g * Isp[0] * md  # Forca propulsiva constante
        elif t <= ts[0]:
            m = m0 - mp[0]  # Entre a queima e a separacao
            ft = 0
        else:
            m = m0 - mp[0] - ms[0]  # Apos a separacao do motor foguete
            ft = 0
        return ft, m

    def propulsor_2_estagios(t, ti, tq, ts, Isp, mp, ms, m0, g):
        if t <= ti[0]:
            m = m0  # Massa inicial
            ft = 0  # Forca propulsiva nula
        elif t <= tq[0]:
            md = -mp[0] / (tq[0] - ti[0])  # Taxa de queima continua
            m = m0 + md * (t - ti[0])  # Esta queimando o primeiro estagio
            ft = -g * Isp[0] * md  # Forca propulsiva constante
        elif t <= ts[0]:
            m = m0 - mp[0]  # Entre a queima e a separacao
            ft = 0
        elif t <= ti[1]:
            m = m0 - mp[0] - ms[0]  # Entre a separacao e a ignicao
            ft = 0
        elif t <= tq[1]:
            md = -mp[1] / (tq[1] - ti[1])  # Taxa de queima continua no segundo estagio
            m02 = m0 - mp[0] - ms[0]
            m = m02 + md * (t - ti[1])  # Durante a queima do segundo estagio
            ft = -g * Isp[1] * md  # Forca propulsiva constante
        elif t <= ts[1]:
            m = m0 - mp[0] - ms[0] - mp[1]  # Apos a queima do segundo estagio e antes da separacao do mesmo
            ft = 0
        else:
            m = m0 - mp[0] - ms[0] - mp[1] - ms[1]  # Apos a separacao do segundo estagio
            ft = 0
        return ft, m

    def propulsor_3_estagios(t, ti, tq, ts, Isp, mp, ms, m0, g):
        if t <= ti[0]:
            m = m0  # Massa inicial
            ft = 0  # Forca propulsiva nula
        elif t <= tq[0]:
            md = -mp[0] / (tq[0] - ti[0])  # Taxa de queima continua
            m = m0 + md * (t - ti[0])  # Esta queimando o primeiro estagio
            ft = -g * Isp[0] * md  # Forca propulsiva constante
        elif t <= ts[0]:
            m = m0 - mp[0]  # Entre a queima e a separacao
            ft = 0
        elif t <= ti[1]:
            m = m0 - mp[0] - ms[0]  # Entre a separacao e a ignicao
            ft = 0
        elif t <= tq[1]:
            md = -mp[1] / (tq[1] - ti[1])  # Taxa de queima continua no segundo estagio
            m02 = m0 - mp[0] - ms[0]
            m = m02 + md * (t - ti[1])  # Durante a queima do segundo estagio
            ft = -g * Isp[1] * md  # Forca propulsiva constante
        elif t <= ts[1]:
            m = m0 - mp[0] - ms[0] - mp[1]  # Apos a queima do segundo estagio e antes da separacao do mesmo
            ft = 0
        elif t <= ti[2]:
            m = m0 - mp[0] - ms[0] - mp[1] - ms[1]  # Entre a separacao e a ignicao
            ft = 0
        elif t <= tq[2]:
            md = -mp[2] / (tq[2] - ti[2])  # Taxa de queima continua no terceiro estagio
            m03 = m0 - mp[0] - ms[0] - mp[1] - ms[1]
            m = m03 + md * (t - ti[2])  # Durante a queima do terceiro estagio
            ft = -g * Isp[2] * md  # Forca propulsiva constante
        elif t <= ts[2]:
            m = m0 - mp[0] - ms[0] - mp[1] - ms[1] - mp[2]  # Apos a queima do terceiro estagio e antes da separacao doDesculpe, a resposta anterior foi cortada. Aqui está a continuação:
            ft = 0
        else:
            # Apos a separacao do terceiro estagio
            m = m0 - mp[0] - ms[0] - mp[1] - ms[1] - mp[2] - ms[2]
            ft = 0
        return ft, m

    def propulsor_3_estagios_2ig(t, ti, tq, ts, Isp, mp, ms, m0, mL, g):
        """
        Cálculo da massa e tração
        """
        # Antes da ignição
        if t <= ti[0]:
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
        elif t <= ti[2]:
            # Entre a separação e a ignição
            m = m0 - mp[0] - ms[0] - mp[1] - ms[1]
            ft = 0
        elif t <= tq[2]:
            # Taxa de queima contínua no terceiro estágio - primeira ignicao
            md = -mp[2] / (tq[2] - ti[2])
            # Durante a queima do terceiro estágio
            m03 = m0 - mp[0] - ms[0] - mp[1] - ms[1]
            m = m03 + md * (t - ti[2])
            # Força propulsiva constante
            ft = -g * Isp[2] * md
        elif t <= ti[3]:
            # Antes da nova queima do terceiro estagio
            m = m0 - mp[0] - ms[0] - mp[1] - ms[1] - mp[2]
            ft = 0
        elif t <= tq[3]:
            # Taxa de queima contínua no terceiro estágio - segunda ignicao
            md = -mp[3] / (tq[3] - ti[3])
            # Durante a queima do terceiro estágio
            m03 = m0 - mp[0] - ms[0] - mp[1] - ms[1] - mp[2]
            m = m03 + md * (t - ti[3])
            # Força propulsiva constante
            ft = -g * Isp[2] * md
        elif t <= ts[2]:
            # Após a queima do terceiro estágio e antes da separação do mesmo
            m = m0 - mp[0] - ms[0] - mp[1] - ms[1] - mp[2] - mp[3]
            ft = 0
        else:
            # Após a separação do terceiro estágio
            m = mL
            ft = 0

        return ft, m
    

    # As variáveis propulsivas dependem do estágio atual
    # Número de estágios
    
    N = len(ent.ti)
    if N == 1:
        ft, m = propulsor_1_estagio(t, ent.ti, ent.tq, ent.ts, ent.Isp, ent.mp, ent.ms, ent.m0, ent.g)
    elif N == 2:
        ft, m = propulsor_2_estagios(t, ent.ti, ent.tq, ent.ts, ent.Isp, ent.mp, ent.ms, ent.m0, ent.g)
    elif N == 3:
        ft, m = propulsor_3_estagios(t, ent.ti, ent.tq, ent.ts, ent.Isp, ent.mp, ent.ms, ent.m0, ent.g)
    else: 
        ft, m = propulsor_3_estagios_2ig(t, ent.ti, ent.tq, ent.ts, ent.Isp, ent.mp, ent.ms, ent.m0, ent.mL, ent.g)    




    # Para altitudes acima de 200km, alinha o vetor de tração com a
    # velocidade inercial ao inves da relativa
    # Desmembra o vetor de estado
    

    # Altitude
    h = r - ent.Requat
    if h < 200e3:
        epsl = 0
        mu = 0  # Tração alinhada com a velocidade relativa
    else:
        # Vetor velocidade inercial
        _, phii, Ai = Vrel2Vine(V, phi, A, ent.we, r, delta)
        
        # Angulos propulsivos para que a tração seja alinhada com a
        # velocidade inercial
        mu = asin(np.cos(A) * np.cos(phii) * np.sin(Ai) - np.sin(A) * np.cos(phii) * np.cos(Ai))
        epsl = -atan2(-np.cos(phi) * np.sin(phii) + np.sin(phi) * np.sin(A) * np.cos(phii) * np.sin(Ai) + \
                      np.sin(phi) * np.cos(A) * np.cos(phii) * np.cos(Ai),
                      np.sin(phi) * np.sin(phii) + \
                      np.cos(phi) * np.sin(A) * np.cos(phii) * np.sin(Ai) + np.cos(phi) * np.cos(A) * np.cos(phii) * np.cos(Ai))

    return ft, m, mu, epsl



