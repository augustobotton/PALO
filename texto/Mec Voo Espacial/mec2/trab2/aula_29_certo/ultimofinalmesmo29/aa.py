import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from scipy.constants import g
from scipy.integrate import solve_ivp

###################################################################################################################################
###################################################################################################################################
###################################################################################################################################

global we, Requat, lc, dT, h0, l_trilho, ti, tq, ts, Isp, mp, ms, m0, mL, sinalPhii, achouApogeu

###################################################################################################################################
###################################################################################################################################
###################################################################################################################################


# As funções "propulsao_N_estagios", "atm_padrao", "aerodinamica_N_estagios", "grav_axisimetrico" 
# e "parametros_manobra_adquire_gso" precisam ser definidas. Os nomes dessas funções são traduções diretas do MATLAB, mas 
# as implementações dessas funções dependem do que elas estão fazendo exatamente no seu código original.

def dinamica_foguete(t, X):
    # global we, Requat, lc, dT, h0, l_trilho, rho
    # Vetor de estado
    V, A, phi, r, delta = X[0], X[1], X[2], X[3], X[4]



    if V < 0:
        V = 0  # Evita velocidade negativa 

    # Função para cálculo da massa e tração em função do tempo
    ft, m, mu, epsl = propulsao_N_estagios(t, X)
    
    # Função para cálculo do modelo atmosférico
    h = r - Requat  # Altitude
    T, _, _, rho, _, Mach, mu, _, Kn, _, _,R = atm_padrao(h, V, lc, dT)
    
    # Cálculo do modelo aerodinâmico
    D, fy, L = aerodinamica_N_estagios(t, V, h, Mach, Kn, T, R,rho)
    
    # Cálculo da gravidade
    gc, gd = grav_axisimetrico(r, delta)

    # Equações de cinemática de translação
    rp = V * np.sin(phi)
    deltap = (V / r) * np.cos(phi) * np.cos(A)
    lonp = (V * np.cos(phi) * np.sin(A)) / (r * np.cos(delta))

    # Equações de dinâmica de translação
    Vp = (1 / m) * (ft * np.cos(epsl) * np.cos(mu) - D - m * gc * np.sin(phi) + m * gd * np.cos(phi) * np.cos(A) - 
        m * we**2 * r * np.cos(delta) * (np.cos(phi) * np.cos(A) * np.sin(delta) - np.sin(phi) * np.cos(delta)))
    Ap = (1 / (m * V * np.cos(phi))) * (m * (V**2 / r) * np.cos(phi)**2 * np.sin(A) * np.tan(delta) + ft * np.sin(mu) + fy 
        - m * gd * np.sin(A) + m * we**2 * r * np.sin(A) * np.sin(delta) * np.cos(delta) - 2 * m * we * V * 
        (np.sin(phi) * np.cos(A) * np.cos(delta) - np.cos(phi) * np.sin(delta)))
    phip = (1 / (m * V)) * (m * (V**2 / r) * np.cos(phi) + ft * np.sin(epsl) * np.cos(mu) + L - m * gc * np.cos(phi) - 
        m * gd * np.sin(phi) * np.cos(A) + m * we**2 * r * np.cos(delta) * (np.sin(phi) * np.cos(A) * np.sin(delta) + 
        np.cos(phi) * np.cos(delta)) + 2 * m * we * V * np.sin(A) * np.cos(delta))

    # Saturação da altitude
    if h < 0:  # Altitude negativa não é permitida
        # Mantém as derivadas nulas
        rp = 0; deltap = 0; lonp = 0; Vp = 0; Ap = 0; phip = 0
    
    # Modela o trilho de lançamento
    H = h - h0  # Altura
    if (H <= l_trilho) and (t <= 10):  # Verifica se a altura é menor que l_trilho nos primeiros segundos da simulação
        Ap = 0; phip = 0  # Anula as derivadas dos ângulos de orientação da velocidade

    # Modela a navegação para o segundo disparo do motor do terceiro estágio
    parametros_manobra_adquire_gso(t, m, X)  # Atualiza variáveis globais

    # Derivada do vetor de estado
    Xp = np.array([Vp, Ap, phip, rp, deltap, lonp])
  
    return Xp


###################################################################################################################################
###################################################################################################################################
###################################################################################################################################

from math import asin, atan2, cos, sin

def propulsao_N_estagios(t, X):
    
    # global ti, tq, ts, Isp, mp, ms, m0, g, mL, we, Requat
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
    else:
        ft, m = propulsor_3_estagios_2ig(t, ti, tq, ts, Isp, mp, ms, m0, mL, g)
       
    
    # Para altitudes acima de 200km, alinha o vetor de tração com a velocidade inercial ao invés da relativa
    # Desmembra o vetor de estado
    V, A, phi, r, delta = X[0], X[1], X[2], X[3], X[4]
    
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

###################################################################################################################################
###################################################################################################################################
###################################################################################################################################


def atm_padrao(h, v, lc, dT):

    # Funcao para calculo da atmosfera padrao para a altitude geometrica de
    # zero ate 2000 km
    # Eh utilizado o modelo apresentado na referencia
    # TEWARI, A. Atmospheric and Space Flight Dynamics:
    # Modelling and simulation with MATLAB and Simulink. Boston: Birkhauser, 2007.
    # Capitulo 9 - Planetary atmosphere
    # O modelo possui 21 camadas, sendo uma combinacao dos modelos de atmosfera
    # padrao norte americanos de 1962 e 1976.
    # Todas as camadas consideradas possuem variacao linear de temperatura ou
    # sao isotermicas. Na faixa de altitude de zero ate 86km, eh utilizado o
    # perfil de temperaturas da atmosfera padrao norte americana do ano de
    # 1976. Para a faixa de altitude de 86km a 2000km, eh utilizado o modelo de 1962.
    # Apesar de a atmosfera nao apresentar equilibrio termico e quimico acima
    # de 86km, sendo os perfis de temperatura, nessas regioes, nao lineares,
    # esta funcao adota do modelo linear de temperaturas de 1962 como uma aproximacao.
    # Entradas:
    # h - altitude geometrica [m]
    # v - velocidade do veiculo em relacao ao escoamento nao perturbado, em
    # metros por segundo [m/s]
    # lc - comprimento caracteristico do veiculo, em metros [m]
    # Saidas
    # T - Temperatura cinetica, em Kelvin [K]
    # Tm - Temperatura de escala molecular, em Kelvin [K]
    # p - Pressao, em Pascal [N/m**2]
    # rho - Densidade [kg/m**2]
    # ainf - Velocidade do som [m/s]
    # M - Numero de Mach [adm]
    # mu - Coeficiente de viscosidade dinamica [kg/m*s]
    # Pr - Numero de Prandtl [adm]
    # Kn - Numero de Knudsen [adm]
    # d - Parametro de regime de escoamento [adm]
    # Re - Numero de Reynolds [adm]
    # Entrada de dados
    # Vetores com os seguintes dados, cada linha eh uma camada do modelo de
    # atmosfera padrao: altitude no inicio da camada (m), temperatura no inicio
    # da camada (K), constante de gas ideal do ar na camada (J/kg.K), taxa de
    # lapso termico (K/m), peso molecular M
    hi = 1e3*np.array([0, 11.0191, 20.0631, 32.1619, 47.3501, 51.4125, 71.8020, 86, 100, 110, 120, 150,
                   160, 170, 190, 230, 300, 400, 500, 600, 700])  # km
      # m
    Ti = np.array([288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65, 186.946, 210.65, 260.65, 360.65,
                   960.65, 1110.65, 1210.65, 1350.65, 1550.65, 1830.65, 2160.65, 2420.65, 2590.65, 2700.65])
    Ri = np.array([287.0, 287.0, 287.0, 287.0, 287.0, 287.0, 287.02, 287.02, 287.84, 291.06, 308.79, 311.80,
                   313.69, 321.57, 336.68, 366.84, 416.88, 463.36, 493.63, 514.08, 514.08])
    a = 1e-3*np.array([-6.5, 0.0, 1.0, 2.8, 0.0, -2.8, -2.0, 1.693, 5.0, 10.0, 20.0, 15.0, 10.0, 7.0, 5.0,
                  4.0, 3.3, 2.6, 1.7, 1.1, 0.0])  # °/km
    
    
    Mi = np.array([28.9644, 28.9644, 28.9644, 28.9644, 28.9644, 28.9644, 28.9644, 28.9644,
                   28.88, 28.56, 28.07, 26.92, 26.66, 26.4, 25.85, 24.7, 22.66, 19.94, 17.94, 16.84, 16.17])

    # Constantes
    M0 = Mi[0]  # Peso molecular ao nivel do mar
    g0 = 9.80665  # Valor ao nivel do mar da aceleracao da gravidade (m/s**2)
    Na = 6.0220978e23  # Numero de Avogadro
    sigma = 3.65e-10  # Diametro de colisao para o ar (m)
    m0 = 28.964e-3  # Massa molar do ar ao nivel do mar (kg/Mol)
    P0 = 1.01325e5  # Pressao padrao ao nivel do mar (N/m**2)
    Requat = 6378.14e3  # Raio medio da Terra (m)
    gamma = 1.405  # Razao de calores especificos ao nivel do mar
    # Constantes calculadas
    # Numero beta associado a distancia radial media do nivel do mar
    beta = 2/Requat



    
    # Identifica a camada a qual a altitude pertence

    # Contador
    i = 0
    if h < 0:
        # Altitude negativa. Os resultados apresentados dizem respeito a h=0.
        i = 0; h = 0
    elif h > 2000e3:
        # A altitude fornecida esta acima do limite superior de 2.000 km.
        # Os resultados apresentados dizem respeito a h=2.000 km.
        i = 20; h = 2000e3
    else:
        for i in range(21):
            if (i == 20):
                break
            elif ((h >= hi[i]) and (h < hi[i+1])):
                break

    # Realiza os calculos
    # Pressao no inicio da camada i
    Pi = P0  # Pressao inicial da camada - inicializa com o valor ao nivel do mar
    px = P0  # Variavel auxiliar para guardar o valor de pressao inicial na camada anterior
    for j in range(1, i+1):
        if a[j-1] != 0:  # Verifica se a camada nao eh isotermica
        # Calcula a pressao inicial da camada j a partir do modelo da
        # camada j-1
            A = 1+(a[j-1]*(hi[j]-hi[j-1]))/Ti[j-1]
            B = -(g0/(Ri[j]*a[j-1]))*(1+beta*((Ti[j-1]/a[j-1])-hi[j-1]))
            C = (g0*beta/(Ri[j]*a[j-1]))*(hi[j]-hi[j-1])
            Pi = px*(A**B)*np.exp(C)
            px = Pi  # O valor atual sera o valor anterior na proxima iteracao
        else:
            # Calcula a pressao inicial da camada i pelo modelo isotermico
           Pi = px*np.exp(-(g0/(Ri[j]*Ti[j-1]))*(hi[j]-hi[j-1])*(1-(beta/2)*(hi[j]-hi[j-1])))
           px = Pi  # O valor atual sera o valor anterior na proxima iteracao
    # Temperatura padrao (molecular) - usada nos calculos internos
    Tm = Ti[i]+a[i]*(h-hi[i])
    Tm=Tm+dT

    # Interpola o valor de R e M
    
    if i >= 20:
        R = Ri[i]
        M = Mi[i]
    else:
        R = Ri[i]+((Ri[i+1]-Ri[i])/(hi[i+1]-hi[i]))*(h-hi[i])
        M = Mi[i]+((Mi[i+1]-Mi[i])/(hi[i+1]-hi[i]))*(h-hi[i])
        # Temperatura padrao (cinetica) - saida da funcao
    T = (M/M0)*Tm
    # Pressao
    if a[i] != 0:  # Verifica se a camada nao eh isotermica
        # Calcula a pressao
        A = 1+(a[i]*(h-hi[i]))/Ti[i]
        B = -(g0/(R*a[i]))*(1+beta*((Ti[i]/a[i])-hi[i]))
        C = (g0*beta/(R*a[i]))*(h-hi[i])
        p = Pi*(A**B)*np.exp(C)
    else:
        # Calcula a pressao pelo modelo isotermico
        p = Pi*np.exp(-(g0/(R*Ti[i]))*(h-hi[i])*(1-(beta/2)*(h-hi[i])))
        
    # Calcula a densidade
    rho = p/(R*Tm)
    # Velocidade do som
    ainf = np.sqrt(gamma*R*Tm)
    # Numero de Mach
    Mach = v/ainf
    # Coeficiente de viscosidade dinamica
    mu = 1.458e-6*(Tm**(3/2))/(Tm+110.4)
    # Numero de Prandtl
    cp = R*gamma/(gamma-1)
    kT = (2.64638e-3*(Tm**(3/2)))/(Tm+245.4*(10**(-12/Tm)))
    Pr = mu*cp/kT
    # Numero de Knudsen
    lam = m0/(np.sqrt(2)*np.pi*sigma**2*rho*Na)
    Kn = lam/lc
    # Parametro de regime de escoamento
    if Kn >= 10:
        d = 1
    elif Kn <= 0.01:
        d = 2
    else:
        d = 3
    # Re - Numero de Reynolds [adm]
    Reynolds = rho*v*lc/mu
    return T, Tm, p, rho, ainf, Mach, mu, Pr, Kn, d, Reynolds,R  



###################################################################################################################################
###################################################################################################################################
###################################################################################################################################

def aerodinamica_N_estagios(t, V, h, Mach, Kn, T, R,rho):

    def area1estagio(t, ts, Sr):
        if t <= ts[0]:
            # Foguete e carga útil
            S = Sr[0]
        else:
            S = Sr[1]  # Carga útil
        return S

    def area2estagios(t, ts, Sr):
        if t <= ts[0]:
            # Todos os estágios
            S = Sr[0]
        elif t <= ts[1]:
            # Segundo estágio e carga útil
            S = Sr[1]
        else:
            # Carga útil
            S = Sr[2]
        return S

    def area3estagios(t, ts, Sr):
        if t <= ts[0]:
            # Todos os estágios
            S = Sr[0]
        elif t <= ts[1]:
            # Segundo estágio
            S = Sr[1]
        elif t <= ts[2]:
            # Terceiro estágio e carga útil
            S = Sr[2]
        else:
            # Carga útil
            S = Sr[3]
        return S

    # global rho
    # Coeficiente de arrasto em função do número de Mach e de Knuden
    CD = modelo_aerodinamico(V, h, Mach, Kn, T, R)
    # Fator de correção do arrasto a partir de dados de túnel de vento
    CD = fc * CD
    
    # A área de referência depende do estágio atual
    # Número de estágios
    N = len(ts)
    if N == 1:
        S = area1estagio(t, ts, Sr)
    elif N == 2:
        S = area2estagios(t, ts, Sr)
    else:
        S = area3estagios(t, ts, Sr)
 
    # Forças
    D = 0.5 * rho * V**2 * S * CD
    fy, L = 0, 0
    
    return D, fy, L

###################################################################################################################################
###################################################################################################################################
###################################################################################################################################


from scipy.interpolate import pchip_interpolate
from math import sqrt, pi, log10, sin
def modelo_aerodinamico(V,h,Mach,Kn,T,R):
    # Data for the continuous drag coefficient calculation
    CD = np.array([[0.0000,    0.4736],
        [0.1609,    0.4736],
        [0.3448,    0.4764],
        [0.4828,    0.4875],
        [0.5261,    0.5001],
        [0.6424,    0.5141],
        [0.7091,    0.5389],
        [0.7521,    0.5646],
        [0.8463,    0.5905],
        [0.8736,    0.6153],
        [0.9065,    0.6667],
        [0.9378,    0.6375],
        [0.9492,    0.7029],
        [0.9530,    0.7217],
        [0.9721,    0.7858],
        [0.9759,    0.8046],
        [0.9949,    0.8680],
        [0.9987,    0.8866],
        [1.0049,    1.0313],
        [1.0115,    0.9833],
        [1.0178,    0.9502],
        [1.0247,    1.0133],
        [1.0345,    0.9667],
        [1.1264,    1.0125],
        [1.1494,    0.9958],
        [1.2184,    0.9833],
        [1.2644,    0.9653],
        [1.2656,    0.9569],
        [1.3037,    0.9398],
        [1.3494,    0.9193],
        [1.3952,    0.8998],
        [1.4409,    0.8797],
        [1.4866,    0.8597],
        [1.5324,    0.8415],
        [1.5862,    0.8194],
        [1.6696,    0.7925],
        [1.7382,    0.7706],
        [1.8068,    0.7491],
        [1.8754,    0.7277],
        [1.9770,    0.6944],
        [2.1842,    0.6597],
        [2.2757,    0.6445],
        [2.3557,    0.6315],
        [2.6073,    0.5941],
        [2.7331,    0.5785],
        [2.8703,    0.5596],
        [2.9885,    0.5458],
        [3.1695,    0.5396],
        [3.4306,    0.5272],
        [3.6822,    0.5171],
        [3.9719,    0.5109],
        [4.1853,    0.5060],
        [4.4369,    0.5019],
        [4.6885,    0.4979],
        [4.9545,    0.4920],
        [5.1821,    0.4923],
        [5.4432,    0.4921],
        [5.6948,    0.4921],
        [5.9463,    0.4908],
        [6.1979,    0.4908],
        [6.4495,    0.4901],
        [6.7011,    0.4894],
        [6.9526,    0.4893],
        [7.2042,    0.4881],
        [7.4558,    0.4880],
        [7.7073,    0.4872],
        [7.9589,    0.4866],
        [8.2105,    0.4865],
        [8.4621,    0.4853],
        [8.7136,    0.4852],
        [8.9652,    0.4845],
        [9.2168,    0.4839],
        [9.4683,    0.4837],
        [9.7199,    0.4825],
        [10.2231,   0.4823],
        [10.4746,   0.4825],
        [10.7262,   0.4825],
        [10.9778,   0.4825],
        [11.2294,   0.4825],
        [11.4809,   0.4825],
        [11.7325,   0.4824]])

    def interpCD(Mach):
        """
        Function to calculate the continuous drag coefficient.
        """
        return pchip_interpolate(CD[:,0], CD[:,1], Mach)


    """
    Drag model according to the reference 
    TEWARI, A. Atmospheric and Space Flight Dynamics: 
    Modelling and simulation with MATLAB and Simulink. Boston: Birkhauser, 2007.
    Example 12.6
    Inputs
    V: speed (m/s)
    M: Mach number
    Kn: Knudsen number
    T: Temperature (K)
    R: Ideal gas constant of air (J/kg.K)
    Output
    CD: Drag coefficient
    """
    s=V/sqrt(2*R*T)
    CDfm=1.75+sqrt(pi)/(2*s)
    CDc=interpCD(Mach)
    if h<2000e3:
        if Kn<0.0146:
            CD=CDc
        elif Kn<14.5:
            CD=CDc+(CDfm-CDc)*((1/3)*log10(Kn/sin(30*pi/180))+0.05113)
        else:
            CD=CDfm
    else:
        CD=0
    return CD

###################################################################################################################################
###################################################################################################################################
###################################################################################################################################

def grav_axisimetrico(r,delta):

    phi=np.pi/2-delta
    
    gr = -(mut) / (r ** 2) - (1 / r) * mut * (-((J2 * Requat ** 2 * (-1 + 3 * np.cos(phi) ** 2)) / (r ** 3)) - (1.5 * J3 * Requat ** 3 * (-3 * np.cos(phi) +
    5 * np.cos(phi) ** 3)) / (r ** 4) -(J4 * Requat ** 4 * (3 - 30 * np.cos(phi) ** 2 + 35 * np.cos(phi) ** 4)) / (2 * r ** 5)) + \
    (mut / (r ** 2)) * ((0.5 * J2 * Requat ** 2 * (-1 + 3 * np.cos(phi) ** 2)) / (r ** 2) + (0.5 * J3 * Requat ** 3 * (-3 * np.cos(phi) + \
    5 * np.cos(phi) ** 3)) / (r ** 3) + (J4 * Requat ** 4 * (3 - 30 * np.cos(phi) ** 2 + 35 * np.cos(phi) ** 4)) / (8 * r ** 4))

    gphi = -(1. / np.square(r)) * mut * (
    -(3. * J2 * np.square(Requat) * np.cos(phi) * np.sin(phi)) / np.power(r, 2)
    + (0.5 * J3 * np.power(Requat, 3) * (3. * np.sin(phi) - 15. * np.square(np.cos(phi)) * np.sin(phi))) / np.power(r, 3)
    + (J4 * np.power(Requat, 4) * (60. * np.cos(phi) * np.sin(phi) - 140. * np.power(np.cos(phi), 3) * np.sin(phi))) / (8. * np.power(r, 4)))

    gc = -gr
    gdel=-gphi

    return gc, gdel
 

###################################################################################################################################
###################################################################################################################################
###################################################################################################################################


def parametros_manobra_adquire_gso(t, m, X):
    global achouApogeu, sinalPhii, Tq32
    # Desmembra o vetor de estado
    
    V = X[0]
    A = X[1]
    phi = X[2]
    r = X[3]
    delta = X[4]
    #lon = X[5]
    
    # Vetor velocidade inercial
    Vi, phii, _ = Vrel2Vine(V, phi, A, we, r, delta)
    
    # Realização de uma sequência de testes para verificar a ocorrência do
    # apogeu da órbita GTO. Quando ele ocorre, determina os parâmetros da manobra
    if r > 0.9 * agso:  # Só inicia a verificação quando está perto do apogeu da órbita GTO    
        if not achouApogeu:  # Se não achou o apogeu, entra na rotina de busca
            if np.sign(phii) != sinalPhii:  # Se o sinal for diferente, phii passou por zero, o foguete chegou no apogeu
                achouApogeu = True  # Encontrou o apogeu, ao setar essa variável, só vai entrar aqui uma vez
                ti[3] = t  # Tempo da segunda ignição do motor do terceiro estágio
                
                # Calcula o tempo de duração da queima para propiciar o DeltaV necessário
                DVgso = vgso - Vi  # Quando se passa nos testes acima, "vi" é a velocidade de apogeu da GTO
                mp32 = (m * np.exp(DVgso / (Isp[2] * g)) - m) / np.exp(DVgso / (Isp[2] * g))  # Massa de propelente necessária
                Tq32 = Tq3 * mp32 / mp3  # Duração da queima necessária
                
                # Verifica se o tempo é maior que o máximo, se ocorrer, corrige
                if Tq31 + Tq32 > Tq3:
                    Tq32 = Tq3 - Tq31
                
                tq[3] = ti[3] + Tq32  # Tempo de fim de queima
                ts[2] = tq[3] + Ts3  # Tempo de separação
    
    sinalPhii = np.sign(phii)  # Guarda o sinal de phi inercial para verificar mudança na próxima iteração
    
    return sinalPhii, achouApogeu, ti, tq, ts, Tq32


###################################################################################################################################
###################################################################################################################################
###################################################################################################################################


def Vrel2Vine(vr, phir, Ar, we, r, dt):
    """
    Função para converter a velocidade, elevação e azimute da velocidade
    relativa para a velocidade inercial

    Entradas:
    vr (m/s): velocidade relativa
    phir (rad): inclinação da velocidade relativa
    Ar (rad): azimute da velocidade relativa
    we (rad/s): velocidade de rotação do referencial girante
    r (m): distância radial até a origem do referencial inercial
    dt (rad): latitude

    Saídas:
    vi (m/s): magnitude da velocidade com respeito ao referencial inercial
    phii (rad): ângulo de elevação da velocidade inercial (angulo de trajetoria)
    Ai (rad): ângulo de azimute da velocidade inercial
    """

    # Cálculos
    Ai = np.arctan2(vr * np.cos(phir) * np.sin(Ar) + we * r * np.cos(dt), vr * np.cos(phir) * np.cos(Ar))
    Ai = Ai if Ai >= 0 else Ai + 2 * np.pi

    vi = np.sqrt(vr*2 + 2 * vr * np.cos(phir) * np.sin(Ar) * r * we * np.cos(dt) + r**2 * we*2 * np.cos(dt)**2)

    phii = np.arctan2(np.sin(phir) * np.cos(Ai), np.cos(phir) * np.cos(Ar))
    if abs(phii) > np.pi/2:
        phii = phii + np.pi if phii < np.pi/2 else phii - np.pi

    return vi, phii, Ai


###################################################################################################################################
###################################################################################################################################
###################################################################################################################################

def long_ECEF2ECI(t, long, we, tg):
    """
    Função para calcular a longitude celeste a partir da longitude
    fixa ao planeta
    
    Parâmetros:
    t (s) - Tempo no qual se deseja saber a longitude celeste
    long (rad) - Longitude relativa ao referencial fixo ao planeta
    we (rad/s) - Velocidade de rotação do planeta
    tg (s) - Tempo no qual o meridiano de referência tem longitude celeste
    nula
    
    Retorno:
    long_c (rad) - Longitude celeste no tempo t
    """
    long_c = long + we * (t - tg)
    return long_c


###################################################################################################################################
###################################################################################################################################
###################################################################################################################################

def RvelPolar2RvelRet(v, A, phi, r, lat, long):
    """
    Função para converter velocidade do sistema LVLH (coordenadas polares)
    para o sistema ECI ou ECEF retangular
    A velocidade pode ser a relativa ou a inercial, o resultado final será
    correlato

    Entradas:
    v (m/s): Módulo do vetor velocidade
    A (rad): Azimute da velocidade
    phi (rad): Elevação da velocidade
    r (m): distância radial
    lat (rad): Latitude 
    long (rad): Longitude no referencial desejado (ECI ou ECEF)

    Saídas:
    R (m): Vetor posição em coordenadas retangulares no
    sistema ECI ou ECEF (dependendo da entrada de dados de longitude)
    V (m/s): Vetor velocidade em coordenadas retangulares no
    sistema ECI ou ECEF (dependendo da entrada de dados de longitude). Pode
    ser a velocidade relativa ou a inercial, dependendo dos dados de
    velocidade fornecidos.
    """
    # Matriz de conversão do sistema ECI ou ECEF para o LVLH 
    CLH = np.array([[np.cos(lat)*np.cos(long), np.cos(lat)*np.sin(long), np.sin(lat)],
                    [-np.sin(long), np.cos(long), 0],
                    [-np.sin(lat)*np.cos(long), -np.sin(lat)*np.sin(long), np.cos(lat)]])

    # Vetor velocidade em coordenadas cartezianas no sistema LVLH
    Vlvlh = v * np.array([np.sin(phi),
                          np.cos(phi)*np.sin(A),
                          np.cos(phi)*np.cos(A)])

    # Transformação da velocidade para o sistema ECI ou ECEF em coordenadas retangulares
    V = np.dot(CLH.T, Vlvlh)

    # Vetor posição no sistema LHVLH
    Rlvlh = np.array([r, 0, 0])

    # Transformação da posição para o sistema ECI ou ECEF em coordenadas retangulares
    R = np.dot(CLH.T, Rlvlh)

    return R,V


###################################################################################################################################
###################################################################################################################################
###################################################################################################################################
import math

def det_orbita(t0, rc0, vc0, mu):
    """
    Função para determinar parâmetros orbitais a partir de uma observação de
    posição e outra de velocidade, sendo as mesmas tomadas em relação ao
    primário em um problema de dois corpos e escritas no referencial
    celestial.
    """
    # Cálculos
    # Distância radial ao primário no instante observado
    r0 = np.linalg.norm(rc0)
    
    # Vetor quantidade de movimento angular específica no referencial celeste
    hc = np.cross(rc0, vc0)
    
    # Vetor excentricidade no sistema celeste
    ec = np.cross(vc0, hc) / mu - rc0 / r0
    
    # Excentricidade da órbita
    e = np.linalg.norm(ec)
    
    # Módulo do vetor hc
    h = np.linalg.norm(hc)
    
    # Parâmetro da órbita dada
    p = h**2 / mu
    
    # Semi eixo maior
    a = p / (1 - e**2)
    
    # Vetor parâmetro no referencial celeste
    pc = p * np.cross(hc, ec) / (h * e)
    
    # Anomalia verdadeira
    costheta = (p - r0) / (e * r0)
    sintheta = np.dot(rc0, pc) / (r0 * p)
    theta = math.atan2(sintheta, costheta)
    
    # O tempo de perigeu depende do tipo de órbita
    if (0 <= e) and (e < 1):
        tipo = 'e'  # órbita elíptica
    elif e == 1:
        tipo = 'p'  # órbita parabólica
    else:
        tipo = 'h'  # órbita hiperbólica
    
    # Tempo de perigeu
    if tipo == 'e':  # órbita elíptica
        # Movimento médio
        n = math.sqrt(mu / a**3)
        # Anomalia excêntrica
        E = 2 * math.atan(math.sqrt((1 - e) / (1 + e)) * math.tan(theta / 2))
        tau = t0 - (E - e * math.sin(E)) / n
    elif tipo == 'p':  # órbita parabólica
        tau = -((math.tan(theta / 2))*3 + 3 * math.tan(theta / 2)) / (mu / p*3)**(1 / 6)
    else:  # órbita hiperbólica
        # Movimento médio hiperbólico
        n = math.sqrt(-mu / a**3)
        # Anomalia hiperbólica
        H = 2 * math.atanh(math.sqrt((e - 1) / (1 + e)) * math.tan(theta / 2))
        tau = -(e * math.sinh(H) - H) / n
    
    # Linha dos nodos
    # Vetor unitário ao longo do vetor h (no sistema celeste)
    ih = hc / h
    # Vetor unitário ao longo da linha dos nodos (no sistema celeste)
    Kc = np.array([0, 0, 1])
    nc = np.cross(Kc, ih) / np.linalg.norm(np.cross(Kc, ih))
    # Ascensão reta do nodo ascendente
    OMEGA = math.atan2(nc[2], nc[1])
    # Inclinação
    i = math.acos(np.dot(ih, Kc))
    # Vetor unitário ao longo do vetor excentricidade (no referencial celeste)
    ie = ec / e
    # Argumento de perigeu
    cosomega = np.dot(ie, nc)
    sinomega = np.dot(ih, np.cross(nc, ie))
    omega = math.atan2(sinomega, cosomega)
    
    # Vetor de parâmetros de saída
    par_orb = np.array([a, e, tau, OMEGA, i, omega])
    return par_orb


###################################################################################################################################
###################################################################################################################################
###################################################################################################################################

Isp = np.array([251, 271, 315])  # s - Impulso específico dos estagios
mp = np.array([5.5262e+04, 11058, 243.6, 0])  # kg - Massa de propelente dos estagios
Tq1 = 62  # s - NAO MUDA - DADO DOS MOTORES DO 1° ESTAGIO
Tq2 = 64.62  # s - NAO MUDA - DADO DOS MOTORES DO 2° ESTAGIO
Tq3 = 301  # s - TEMPO DE QUEIMA DO 3° ESTAGIO SE ELE IGNITASSE SO UMA VEZ

# Parâmetros de massa estrutural e de carga util
ms = np.array([7750, 1367, 64.7544])  # kg - Massa estrutural dos estagios
mL = 13  # kg - Massa da carga util

# Parâmetros aerodinâmicos e ambientais
fc = 1.28  # Fator de correção do arrasto a partir de dados de tunel de vento
S1 = 4.6*5/3  # m^2 - Area aproximada da seção transversal do primeiro estagio
S2 = 1.5  # m^2 - Area aproximada da seção longitudinal do segundo estagio
S3 = 1.5  # m^2 - Area aproximada da seção longitudinal do terceiro estagio
SL = 1.5  # m^2 -  Area aproximada da seção longitudinal da carga util
lt = 7.33 + 7.1 + 6.28  # m - Comprimento total
l2 = 7.1 + 6.28  #  Comprimento sem o primeiro estagio
l3 = 6.28  #  Comprimento sem o segundo estagio
l4 = 1  # Comprimento da carga util
f2 = (l2/lt)*0.5 + 0.5  # Fator de correção do segundo estagio
f3 = (l3/lt)*0.5 + 0.5  # Fator de correção do terceiro estagio
f4 = (l4/lt)*0.5 + 0.5  # Fator de correção da carga util
Sr = np.array([S1, S2*f2, S3*f3, SL*f4])  # Vetor de areas de referencia para calculo do arrasto
lc = 1.5  # Comprimento característico - diametro dos estagios 2 e superiores
dT = 10  # K - Delta T em relação à atmosfera padrão (que é 15°C no nível do mar)

# Parâmetros da Terra - modelo axis simetrico (WGS-84)
Requat = 6378.1370e3  # m - Raio equatorial da Terra
we = 7.2921150e-5  # (rad/s) - Velocidade inercial de rotação da Terra
g = 9.80665  # m/s^2 - aceleracao da gravidade padrao ao nivel do mar
mut = 3.986004418e14  # m3.s^-2
J2 = 0.00108263  # Constante de Jeffery
J3 = -0.00000254  # Constante de Jeffery
J4 = -0.00000161  # Constante de Jeffery
tg = 0  # s - Tempo em que o meridiano de referência tem longitude celeste nula

# Condições iniciais - Centro espacial de Alcântara (do Google Maps)
h0 = 0  # m - Altitude da base de lançamento
delta0 = -2.3267844 * np.pi / 180  # rad - Latitude inicial
lon0 = -44.4111042 * np.pi / 180  # rad - Longitude inicial
l_trilho = lt  # m - igual ao comprimento total do foguete

# Parâmetros da órbita desejada
ingso = 5 * np.pi / 180  # Inclinação
agso = 42.164140e6  # m
vgso = np.sqrt(mut / agso)

# PARAMETROS PROPULSIVOS E TEMPORAIS A DETERMINAR 
# SAO DEFINIDOS PARA PROPICIAR A INSERCAO ORBITAL
Ts1 = 2  # s
Ts2 = 2  # s
Ts3 = 2  # s
TEq2 = 5  # s
TEq3 = 1360 # s
Tq31 = 262
Tq32 = 39
mp31 = mp[2] * Tq31 / Tq3
mp32 = mp[2] * Tq32 / Tq3
mp3 = mp[2]
ti = np.zeros(5)
tq = np.zeros(5)
ts = np.zeros(5)
ti[0] = 0  # s - Tempo da ignicao do primeiro estagio
tq[0] = ti[0] + Tq1  # s - Tempo do fim da queima do est�gio 1
ts[0] = tq[0] + Ts1  # Tempo da separacao do primeiro estagio
ti[1] = ts[0] + TEq2  # s - Tempo da ignicao do segundo estagio
tq[1] = ti[1] + Tq2  # s - Tempo do fim da queima do est�gio 2
ts[1] = tq[1] + Ts2  # Tempo da separacao do segundo estagio
ti[2] = ts[1] + TEq3  # s - Tempo da primeira ignicao do terceiro estagio
tq[2] = ti[2] + Tq31  # s - Tempo do fim da primeira queima do est�gio 3
ti[3] = 1e10  # s - Tempo da segunda ignicao do terceiro estagio
tq[3] = ti[3] + Tq32  # s - Tempo do fim da segunda queima do est�gio 3
ts[2] = tq[3] + Ts3  # Tempo da separacao do terceiro estagio
mp[2] =  mp31 
mp[3] =  mp32 # Complementa o vetor de massas de propelente

# Parametros globais para procurar o apogeu da orbita de transferencia
sinalPhii = 0
achouApogeu = 0

# Parametros calculados
m0 = np.sum(mp) + np.sum(ms) + mL  # Massa inicial do foguete
r0 = Requat + h0  # Distancia radial inicial

# Estudo simplificado pela equacao de foguete
mpx = (mp[0],mp[1], mp3)  # Razão estrutural do primeiro e segundo estágios
print(mpx)
sigma = ms / (ms + mpx)
m01 = m0  # Massa total na decolagem
m02 = ms[1] + mpx[1] + ms[2] + mpx[2] + mL  # Massa total na ignição do segundo estágio
m03 = ms[2] + mpx[2] + mL  # Massa total na ignição do terceiro estágio
lamb = np.zeros(3)
lamb[0] = m02 / m01  # Razão de carga útil do primeiro estágio
lamb[1] = m03 / m02  # Razão de carga útil do segundo estágio
lamb[2] = mL / m03  # Razão de carga útil do terceiro estágio
lambL = np.prod(lamb)  # Razão de carga útil total
ve = g * Isp  # Velocidade de exaustão
Dv = -np.sum(ve * np.log(sigma + (1 - sigma) * lamb))  # Delta v ideal da configuração original

# Mostra dados na tela
print(f'Area de referencia do foguete com primeiro estagio (m^2): {Sr[0]}')
print(f'Area de referencia do foguete com segundo estagio (m^2): {Sr[1]}')
print(f'Area de referencia do foguete com terceiro estagio (m^2): {Sr[2]}')
print(f'Area de referencia da carga util (m^2): {Sr[3]}')
print(f'Massa inicial antes da queima do primeiro estagio - kg: {m01}')
print(f'Massa inicial antes da queima do segundo estagio - kg: {m02}')
print(f'Massa inicial antes da queima do terceiro estagio - kg: {m03}')
print(f'Massa da carga util - kg: {mL}')
print(f'Razoes estruturais: {sigma}')
print(f'Razoes de carga util: {lamb}')
print(f'Velocidades de exaustao - m/s: {ve}')
print(f'Razao de carga util total: {lambL}')
print(f'Impulso de velocidade total ideal - m/s: {Dv}')

Tempo = 6000 #float(input('Informe o tempo da simulacao (s): '))
V0 = 1 #float(input('Informe o valor inicial da velocidade relativa (m/s): '))
phi0 = 85 #float(input('Informe a condicao inicial do angulo de elevacao (graus): '))
phi0 = phi0 * np.pi / 180
y_t = np.cos(ingso) / np.cos(delta0)  # Faz um teste de factibilidade usando a inclinacao da orbita e a latitude inicial
if abs(y_t) > 1:
    print('Nao eh possivel atingir a inclinacao a partir da latitude inicial. Calculando a menor possivel')
    y = np.sign(y_t)
    
Ai_f = np.arcsin(y_t)  # Condicao final do azimute de velocidade inercial

# Estimativa da condicao inicial de azimute de velocidade relativa
rpgto = Requat + 250e3  # Apogeu de uma orbita de transferencia com 250 km de altitude
agto = (agso + rpgto) / 2  # Semi eixo maior de uma orbita de transferencia com 250 km de altitude
vigto = np.sqrt(mut * (2 / rpgto - 1 / agto))  # Velocidade de uma orbita de transferencia com 250km de altitude
A0 = np.arctan(np.tan(Ai_f) - (rpgto * we * np.cos(delta0)) / (vigto * np.cos(Ai_f)))
  
print(f'Condicao final de azimute de velocidade inercial (grau): {Ai_f * 180 / np.pi}')
print(f'Condicao inicial de azimute de velocidade relativa (grau): {A0 * 180 / np.pi}')

# Simulação
# Condição inicial

X0 = [V0, A0, phi0, r0, delta0, lon0]
# Parâmetros para a função de integração
# Simulação
t0 = 0
options = {'rtol': 1e-8, 'atol': 1e-10, 'max_step': 0.5}
resposta_sim = solve_ivp(dinamica_foguete,(t0,Tempo),y0=X0,**options)
t=resposta_sim.t
X = resposta_sim.y

###################################################################################################################################
###################################################################################################################################
###################################################################################################################################

# Pós-processamento
# Parâmetros da órbita GSO requerida
# Velocidade orbital - passada como variável global para a função de dinâmica,
# que deve  calcular o impulso de velocidade de circularização 
print('*** Orbita GSO requerida ***')
print('Raio da orbita GSO (km)', agso/1e3)
print('Velocidade da orbita GSO (km/s)', vgso/1e3)

# Cálculo de outras variáveis
# Número de instantes de tempo
# Magnitude, azimute e elevação da velocidade relativa

N = len(t)
V = resposta_sim.y[0]
A = resposta_sim.y[1]
phi = resposta_sim.y[2]
# Altitude, latitude e longitude no referencial fixo ao planeta

r = resposta_sim.y[3]
h = np.zeros(N)
delta = resposta_sim.y[4]
lon = resposta_sim.y[5]
m = np.zeros(N)  # Massa
ft = np.zeros(N)  # Força propulsiva
mu = np.zeros(N)
epsl = np.zeros(N)  # Ângulos propulsivos
D = np.zeros(N)  # Força de arrasto
q = np.zeros(N)  # Pressão dinâmica
Mach = np.zeros(N)  # Número de Mach
T = np.zeros(N)  # Temperatura
rho = np.zeros(N)  # Densidade
Vi = np.zeros(N)  # Magnitude da velocidade inercial
phii = np.zeros(N)  # Elevação da velocidade inercial
Ai = np.zeros(N)  # Azimute da velocidade inercial
longc = np.zeros(N)  # Longitude celeste
ee = np.zeros(N)  # Energia específica
a = np.zeros(N)  # Semi eixo maior da órbita
e = np.zeros(N)  # Excentricidade da órbita
tau = np.zeros(N)  # Tempo de perigeu
OM = np.zeros(N)  # Ascenção reta do nodo ascendente
in_ = np.zeros(N)  # Inclinação da órbita
om = np.zeros(N)  # Argumento de perigeu
R0 = np.zeros((N, 3))  # Posição no referencial ECI


#import propulsao_N_estagios as pn
#import aerodinamica_N_estagios as an


for i in range(N):
    
    vetor_parametros = V[i], A[i], phi[i], r[i], delta[i]

    # Posicao no referencial PCPF
    h[i] = r[i]- Requat
    # Forca propulsiva, massa e angulos
    ft[i], m[i], mu[i], epsl[i] = propulsao_N_estagios(t[i], vetor_parametros)
    # Parametros atmosfericos
    T[i], _, _, rho[i], _, Mach[i], _, _, Kn, _, _, R = atm_padrao(h[i], V[i], lc, dT)
    # Forcas aerodinamicas
    
    D[i], _, _ = aerodinamica_N_estagios(t[i], V[i], h[i], (Mach[i]), Kn,( T[i]), R,rho[i])
    # Pressao dinamica
    q[i] = 0.5 * rho[i] * V[i]**2
    # Coordenadas da velocidade inercial no referencial LVLH
    Vi[i], phii[i], Ai[i] = Vrel2Vine(V[i], phi[i], A[i], we, r[i], delta[i])
    # Longitude celeste
    longc[i] = long_ECEF2ECI(t[i], lon[i], we, tg)
    # Energia especifica da orbita
    ee[i] = Vi[i]**2 / 2 - mut / r[i]
    # Posicao e velocidade inercial no referencial ICP
    rc0, vc0 = RvelPolar2RvelRet(Vi[i], Ai[i], phii[i], r[i], delta[i], longc[i])
    R0[i, :] = rc0.T
    # Elementos orbitais
    par_orb = det_orbita(t[i], rc0, vc0, mut)
    a[i], e[i], tau[i], OM[i], in_[i], om[i] = par_orb[0], par_orb[1], par_orb[2], par_orb[3], par_orb[4], par_orb[5]

# Analise de orbita
# Altitude e velocidade inercial no fim da queima do terceiro estagio
for i in range(N):
    if t[i] > tq[2]:  
        break

ifq = i - 2  # Python's index is 0-based
# Tempo do fim da queima do terceiro estagio
tfq = t[ifq]
# Velocidade inercial no fim da queima do terceiro estagio
Vfq = np.full(N, Vi[ifq])
# Altitude no fim da queima do terceiro estagio
hfq = np.full(N, h[ifq])
# Periodo da orbita obtida
P = 2 * np.pi * np.sqrt((Requat + hfq[0]) ** 3 / mut)
print('*** Parametros da Orbita Obtida ***')
print('Velocidade no momento da insercao orbital (km/s)', Vfq[0] / 1e3)
print('Altitude no momento da insercao orbital (km)', hfq[0] / 1e3)
print('Distancia radial no momento da insercao orbital (km)', (hfq[0] + Requat) / 1e3)
print('Semi eixo maior (km)', a[ifq] / 1e3)
print('Periodo(min): ', P / 60)
# Raio do perigeu
rp = a[ifq] * (1 - e[ifq])
# Raio do apogeu
ra = a[ifq] * (1 + e[ifq])
print('Raio do perigeu (km): ', rp / 1e3)
print('Raio do apogeu (km): ', ra / 1e3)
print('Altitude do perigeu (km): ', (rp - Requat) / 1e3)
print('Altitude do apogeu (km): ', (ra - Requat) / 1e3)

# Orbita de transferencia geossincrona (GTO) desejada
print('*** Parametros da Orbita GTO requerida ***')
print('Perigeu da orbita GTO requerida (km)')
rpgto = rp
print(rpgto / 1e3)
print('Apogeu da orbita GTO requerida (km)')
ragto = agso
print(ragto / 1e3)
print('Semi eixo maior da orbita GTO requerida (km)')
agto = (ragto + rpgto) / 2
print(agto / 1e3)
print('Velocidade de perigeu da orbita GTO requerida (km/s)')
vpgto = np.sqrt(mut * (2 / rpgto - 1 / agto))
print(vpgto / 1e3)
print('Velocidade de apogeu da orbita GTO requerida (km/s)')
vagto = np.sqrt(mut * (2 / ragto - 1 / agto))
print(vagto / 1e3)

# Geracao de vetores para tracar grafico
ar = np.full(N, agto)  # Semi eixo maior da orbita GTO requerida
Vir = np.full(N, vpgto)  # Velocidade de perigeu da orbita GTO requerida
eer = -mut / (2 * ar)  # Energia especifica da orbita GTO requerida
eegso = np.full(N, -mut / (2 * agso))  # Energia especifica da orbita GSO requerida

# Tempos de operacao do propulsor do terceiro estagio
print('Tempo de espera para disparo do propulsor do 3º estagio apos a separacao do 2º (s)', TEq3)
print('Duracao do primeiro disparo do motor do 3º estagio (s)', Tq31)
print('Duracao do segundo disparo do motor do 3º estagio (s)', Tq32)
print('Momento do segundo disparo do motor do 3ºestagio (s)', ti[3])  # Indexing is 0-based in python
print('Impulso de velocidade requerido para circularizacao da orbita (km/s)')
DVgso = vgso - vagto
print(DVgso / 1e3)
print('Massa de propelente requerida para circularizacao da orbita (kg)')
# Massa de propelente necessaria
mp32 = (m[ifq] * np.exp(DVgso / (Isp[2] * g)) - m[ifq]) / np.exp(DVgso / (Isp[2] * g))
print(mp32)
print('Massa de propelente disponivel para o 3º disparo (kg)', mp3 - mp31)
print('****** PARAMETROS DA ORBITA FINAL ******')
print('Periodo (min)')
P = 2 * np.pi * np.sqrt(a[-1] ** 3 / mut)
print(P / 60)
print('Semi eixo maior (km)', a[-1])
print('Excentricidade', e[-1])
print('Inclinacao (º)', in_[-1] * 180 / np.pi)

# Figure 1
plt.figure()
plt.subplot(231)
plt.plot(t, V, linewidth=2)
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('V (m/s)')

plt.subplot(232)
plt.plot(t, A * 180 / np.pi, linewidth=2)
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('A ( )')

plt.subplot(233)
plt.plot(t, phi * 180 / np.pi, linewidth=2)
plt.plot(tfq, phi[ifq] * 180 / np.pi, '')
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('\phi ( )')

plt.subplot(234)
plt.plot(t, h / 1e3, linewidth=2)
plt.plot(t, hfq / 1e3, '--')
plt.plot(tfq, hfq[0] / 1e3, '*')
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('h (km)')
plt.legend(['altitude', 'altitude no fim da queima do 3º  estagio'])

plt.subplot(235)
plt.plot(t, delta * 180 / np.pi, linewidth=2)
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('\delta ( )')

plt.subplot(236)
plt.plot(t, lon * 180 / np.pi, linewidth=2)
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('l( )')

# Figure 2
plt.figure()
plt.subplot(221)
plt.plot(t, Vi, linewidth=2)
plt.plot(t, Vir, '--')
plt.plot(t, Vfq, '-.')
plt.plot(tfq, Vfq[0], '*')
plt.grid(True)
plt.xlabel('t (s)')
plt.ylabel('V_i (m/s)')
plt.legend(['Velocidade inercial', 'Velocidade de perigeu da orbita GTO requerida', 'Velocidade no fim da queima do terceiro estagio'])

plt.subplot(222)
plt.plot(t, Ai * 180 / np.pi, linewidth=2)
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('A_i ( )')

plt.subplot(223)
plt.plot(t, phii * 180 / np.pi, linewidth=2)
plt.plot(tfq, phii[ifq] * 180 / np.pi, '')
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('\phi_i ( )')

plt.subplot(224)
plt.plot(t, longc * 180 / np.pi, linewidth=2)
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('\lambda ( )')

# Figure 3
plt.figure()
plt.subplot(221)
plt.plot(t, ft, linewidth=2)
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('f_t (N)')

plt.subplot(222)
plt.plot(t, m, linewidth=2)
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('m (kg)')

plt.subplot(223)
plt.plot(t, mu * 180 / np.pi, linewidth=2)
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('\mu ( )')

plt.subplot(224)
plt.plot(t, epsl * 180 / np.pi, linewidth=2)
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('\epsilon ( )')

# Figure 4
plt.figure()
plt.subplot(311)
plt.plot(t, D, linewidth=2)
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('D (N)')

plt.subplot(323)
plt.plot(t, q, linewidth=2)
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('q (N/m^2)')

plt.subplot(324)
plt.plot(t, Mach, linewidth=2)
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('M (-)')

plt.subplot(325)
plt.plot(t, T - 273.15, linewidth=2)
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('T ( C)')

plt.subplot(326)
plt.plot(t, rho, linewidth=2)
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('\rho (kg/m^3)')

# Figure 5
plt.figure()
plt.subplot(311)
plt.plot(t, ee, linewidth=2)
plt.plot(t, eer, '--', linewidth=2)
plt.plot(t, eegso, '--', linewidth=2)
plt.grid(True)
plt.xlabel('t (s)')
plt.ylabel('\epsilon (J/kg)')
plt.legend(['Energia especifica', 'Energia especifica da orbita GTO requerida', 'Energia especifica da orbita GSO requerida'])

plt.subplot(334)
plt.plot(t, a / 1e3, linewidth=2)
plt.plot(t, ar / 1e3, '--')
plt.plot(t, np.full((N,), Requat) / 1e3, '-.')
plt.grid(True)
plt.xlabel('t (s)')
plt.ylabel('a (km)')
plt.legend(['Semi eixo maior', 'Semi eixo maior da orbita GTO requerida', 'Raio da Terra'])

plt.subplot(335)
plt.plot(t, e, linewidth=2)
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('e (-)')

plt.subplot(336)
plt.plot(t, tau, linewidth=2)
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('\tau (s)')

plt.subplot(337)
plt.plot(t, OM * 180 / np.pi, linewidth=2)
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('\Omega ( )')

plt.subplot(338)
plt.plot(t, in_ * 180 / np.pi, linewidth=2)  # renamed 'in' to 'in_' as 'in' is a reserved word in Python
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('i ( )')

plt.subplot(339)
plt.plot(t, om * 180 / np.pi, linewidth=2)
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('\omega (grau)')


import numpy as np
from mayavi import mlab

import numpy as np
from mayavi import mlab

def desenhar_orbita(longitude, latitude, altitude, texture_path):
    # Converter latitude, longitude e altitude para coordenadas cartesianas
    earth_radius = 6371  # Raio médio da Terra em quilômetros
    x = (earth_radius + altitude) * np.cos(latitude) * np.cos(longitude)
    y = (earth_radius + altitude) * np.cos(latitude) * np.sin(longitude)
    z = (earth_radius + altitude) * np.sin(latitude)

    # Adicionar uma dimensão extra aos vetores
    x = np.expand_dims(x, axis=0)
    y = np.expand_dims(y, axis=0)
    z = np.expand_dims(z, axis=0)

    # Criar uma figura
    fig = mlab.figure(size=(800, 600), bgcolor=(0.1, 0.1, 0.1))

    # Adicionar a textura da Terra ao globo
    surface = mlab.mesh(x, y, z, colormap='gist_earth', scalars=z, figure=fig)
    surface.actor.actor.texture = texture_path

    # Adicionar a trajetória do veículo espacial
    mlab.plot3d(x, y, z, color=(1, 0, 0), tube_radius=0.1, figure=fig)

    # Configurar a visualização do globo
    mlab.view(azimuth=0, elevation=90, distance='auto', focalpoint=(0, 0, 0))

    # Mostrar a figura
    mlab.show()

# Vetores de longitude, latitude e altitude em radianos
longitude = lon  # Vetor de longitudes em radianos
latitude = delta  # Vetor de latitudes em radianos
altitude = h  # Vetor de altitudes
texture_path = 'G:\\.shortcut-targets-by-id\\1-35xUduX2YD_n1NJU0ku_vtiEkY6R714\\trab2\\aula_29_certo\\ultimofinalmesmo29\\earth.jpg'  # Caminho para a imagem de textura da Terra

desenhar_orbita(longitude, latitude, altitude, texture_path)


