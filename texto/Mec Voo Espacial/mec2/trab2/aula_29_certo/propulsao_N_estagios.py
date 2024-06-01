from math import asin, atan2, cos, sin
from parametrosglobais import *
from Vrel2Vine import Vrel2Vine

def propulsao_N_estagios(t, X):
    global ti, tq, ts, Isp, mp, ms, m0, g, mL, we, Requat
    """
    Função para cálculo dos parâmetros propulsivos em função do tempo
    Veículo de até 3 estágios com dupla ignição do terceiro estágio
    Entrada:
        t (s): tempo
        X: vetor de estado
        V = X[0] (m/s): módulo do vetor velocidade relativa com respeito ao planeta girante
        A = X[1] (rad): ângulo de azimute do vetor velocidade relativa com respeito ao eixo z (que aponta para o norte) do sistema uen.
        phi = X[2] (rad): ângulo de elevação do vetor velocidade relativa com respeito ao horizonte local (plano yz do referencial uen)
        r = X[3] (m): distância radial até o centro do planeta
        delta = X[4] (rad): latitude com respeito ao plano equatorial do planeta
        lon = X[5] (rad): longitude planetária
    Saídas:
        ft (N): força propulsiva
        m (kg): massa do foguete em função do tempo
        mu, epsl (rad): Ângulos de apontamento da tubeira
    Hipóteses:
        - Em cada estágio, é assumida taxa de queima contínua, ou seja, não há controle da queima e a mesma é assumida uniforme do início ao fim do propelente
        - A tração de cada estágio é assumida como um pulso retangular perfeito, ou seja, quando acionado, o propulsor vai de tração zero até a máxima, permanecendo nesse patamar constante. Ao fim da queima, a tração cai instantaneamente a zero
    """
    # Número de estágios
    N = len(ti)
    
    if N == 1:
        ft, m = propulsor_1_estagio(t, ti, tq, ts, Isp, mp, ms, m0, g)
    elif N == 2:
        ft, m = propulsor_2_estagios(t, ti, tq, ts, Isp, mp, ms, m0, g)
    elif N == 3:
        ft, m = propulsor_3_estagios(t, ti, tq, ts, Isp, mp, ms, m0, g)
    elif N == 4:
        ft, m = propulsor_3_estagios_2ig(t, ti, tq, ts, Isp, mp, ms, m0, mL, g)
    else:
        raise ValueError("Número de estágios inválido")
    
    # Para altitudes acima de 200km, alinha o vetor de tração com a velocidade inercial ao invés da relativa
    # Desmembra o vetor de estado
    V, A, phi, r, delta = X
    
    # Altitude
    h = r - Requat
    
    if h < 200e3:
        epsl = 0
        mu = 0 # Tração alinhada com a velocidade relativa
    else:
        # Vetor velocidade inercial
        _, phii, Ai = Vrel2Vine(V, phi, A, we, r, delta)
        
        # Ângulos propulsivos para que a tração seja alinhada com a velocidade inercial
        mu = asin(cos(A) * cos(phii) * sin(Ai) - sin(A) * cos(phii) * cos(Ai))
        epsl = -atan2(-cos(phi) * sin(phii) + sin(phi) * sin(A) * cos(phii) * sin(Ai) + sin(phi) * cos(A) * cos(phii) * cos(Ai), 
                       sin(phi) * sin(phii) + cos(phi) * sin(A) * cos(phii) * sin(Ai) + cos(phi) * cos(A) * cos(phii) * cos(Ai))
    return ft, m, mu, epsl

# Modelo do propulsor de foguete de 1 estágio
def propulsor_1_estagio(t,ti,tq,ts,Isp,mp,ms,m0,g):
    # Cálculo da massa e tração em foguete com 1 estágio e carga util
    if t<=ti[0]:
        # Antes da ignição
        m=m0  # Massa inicial
        ft=0  # Força propulsiva nula
    elif t<=tq[0]:
        # Taxa de queima contínua
        md=-mp[0]/(tq[0]-ti[0])
        # Está queimando o primeiro estágio
        m=m0+md*(t-ti[0])
        # Força propulsiva constante
        ft=-g*Isp[0]*md
    elif t<=ts[0]:
        # Entre a queima e a separação
        m=m0-mp[0]
        ft=0
    else:
        # Após a separação do motor foguete
        m=m0-mp[0]-ms[0]
        ft=0
    return ft, m

# Modelo dos propulsores de foguete de 2 estágios
def propulsor_2_estagios(t,ti,tq,ts,Isp,mp,ms,m0,g):
    # Cálculo da massa e tração em foguete com 2 estágios e carga util
    if t<=ti[0]:
        # Antes da ignição
        m=m0   # Massa inicial
        ft=0   # Força propulsiva nula
    elif t<=tq[0]:
        # Taxa de queima contínua
        md=-mp[0]/(tq[0]-ti[0])
        # Está queimando o primeiro estágio
        m=m0+md*(t-ti[0])
        # Força propulsiva constante
        ft=-g*Isp[0]*md
    elif t<=ts[0]:
        # Entre a queima e a separação
        m=m0-mp[0]
        ft=0
    elif t<=ti[1]:
        # Entre a separação e a ignição
        m=m0-mp[0]-ms[0]
        ft=0
    elif t<=tq[1]:
        # Taxa de queima contínua no segundo estágio
        md=-mp[1]/(tq[1]-ti[1])
        # Durante a queima do segundo estágio
        m02=m0-mp[0]-ms[0]
        m=m02+md*(t-ti[1])
        # Força propulsiva constante
        ft=-g*Isp[1]*md
    elif t<=ts[1]:
        # Após a queima do segundo estágio e antes da separação do mesmo
        m=m0-mp[0]-ms[0]-mp[1]
        ft=0
    else:
        # Após a separação do segundo estágio
        m=m0-mp[0]-ms[0]-mp[1]-ms[1]
        ft=0
    return ft, m

# Modelo dos propulsores de foguete de 3 estágios
def propulsor_3_estagios(t,ti,tq,ts,Isp,mp,ms,m0,g):
    # Cálculo da massa e tração
    if t<=ti[0]:
        # Antes da ignição
        m=m0   # Massa inicial
        ft=0   # Força propulsiva nula
    elif t<=tq[0]:
        # Taxa de queima contínua
        md=-mp[0]/(tq[0]-ti[0])
        # Está queimando o primeiro estágio
        m=m0+md*(t-ti[0])
        # Força propulsiva constante
        ft=-g*Isp[0]*md
    elif t<=ts[0]:
        # Entre a queima e a separação
        m=m0-mp[0]
        ft=0
    elif t<=ti[1]:
        # Entre a separação e a ignição
        m=m0-mp[0]-ms[0]
        ft=0
    elif t<=tq[1]:
        # Taxa de queima contínua no segundo estágio
        md=-mp[1]/(tq[1]-ti[1])
        # Durante a queima do segundo estágio
        m02=m0-mp[0]-ms[0]
        m=m02+md*(t-ti[1])
        # Força propulsiva constante
        ft=-g*Isp[1]*md
    elif t<=ts[1]:
        # Após a queima do segundo estágio e antes da separação do mesmo
        m=m0-mp[0]-ms[0]-mp[1]
        ft=0
    elif t<=ti[2]:
        # Entre a separação e a ignição
        m=m0-mp[0]-ms[0]-mp[1]-ms[1]
        ft=0
    elif t<=tq[2]:
        # Taxa de queima contínua no terceiro estágio
        md=-mp[2]/(tq[2]-ti[2])
        # Durante a queima do terceiro estágio
        m03=m0-mp[0]-ms[0]-mp[1]-ms[1]
        m=m03+md*(t-ti[2])
        # Força propulsiva constante
        ft=-g*Isp[2]*md
    elif t<=ts[2]:
        # Após a queima do terceiro estágio e antes da separação do mesmo
        m=m0-mp[0]-ms[0]-mp[1]-ms[1]-mp[2]
        ft=0  
    else:
        # Após a separação do terceiro estágio
        m=m0-mp[0]-ms[0]-mp[1]-ms[1]-mp[2]-ms[2]
        ft=0
    return  ft, m

def propulsor_3_estagios_2ig(t, ti, tq, ts, Isp, mp, ms, m0, mL, g):
    # Cálculo da massa e tração
    if t <= ti[0]:
        # Antes da ignição
        m = m0   # Massa inicial
        ft = 0   # Força propulsiva nula
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
        # Taxa de queima contínua no terceiro estágio - primeira ignição
        md = -mp[2] / (tq[2] - ti[2])
        # Durante a queima do terceiro estágio
        m03 = m0 - mp[0] - ms[0] - mp[1] - ms[1]
        m = m03 + md * (t - ti[2])
        # Força propulsiva constante
        ft = -g * Isp[2] * md
    elif t <= ti[3]:
        # Antes da nova queima do terceiro estágio
        m = m0 - mp[0] - ms[0] - mp[1] - ms[1] - mp[2]
        ft = 0
    elif t <= tq[3]:
        # Taxa de queima contínua no terceiro estágio - segunda ignição
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