import numpy as np



def propulsao_N_estagios(t,ti, ts, tq, Isp, mp, ms, m0, g):

    
    
    def propulsor_1_estagio(t):
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

    def propulsor_2_estagios(t):
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

    def propulsor_3_estagios(t):
        
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
            ft = -g * Isp[2] * md  # Força propulsiva constante
        elif t <= ts[2]:
            m = m0 - mp[0] - ms[0] - mp[1] - ms[1] - mp[2]  
            ft = 0
        else:
            # Apos a separacao do terceiro estagio
            m = m0 - mp[0] - ms[0] - mp[1] - ms[1] - mp[2] - ms[2]
            ft = 0
        return ft, m




    # As variaveis propulsivas dependem do estagio atual
    N = len(ts)  # Numero de estagios
    if N == 1:
        ft, m = propulsor_1_estagio(t)
    elif N == 2:
        ft, m = propulsor_2_estagios(t)
    elif N == 3:
        ft, m = propulsor_3_estagios(t)
    
    return ft, m