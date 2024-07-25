import numpy as np
from matplotlib import pyplot as plt

from src.domain.construtorderesultados.plots_orbitais import plota_orbita
# Supondo que você tenha importado as funções de conversão corretamente
from src.domain.modelos.Simulacao import Simulacao
from src.domain.modelos.orbitas.utilidades import calculos_orbitais
from src.domain.modelos.orbitas.utilidades.funcoes_conversao import componentes_vel_relativa_para_inercial, \
    rvel_polar_para_rvel_retangular, converte_longitude_fixa_planeta_para_longitude_celeste


def plotaresultados(resposta_sim, simulacao: Simulacao):
    t, X = resposta_sim
    N = len(t)
    V = X[0]
    A = X[1]
    phi = X[2]
    r = X[3]
    h = np.zeros(N)
    delta = X[4]
    long = X[5]
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
    we = simulacao.planeta.velocidade_inercial_de_rotacao
    mut = simulacao.planeta.mut
    Requat = simulacao.planeta.raio_equatorial
    agso = simulacao.orbita_alvo.semi_eixo_maior
    vgso = np.sqrt(mut / agso)

    for i in range(N):
        vetor_parametros = V[i], A[i], phi[i], r[i], delta[i], long[i]
        # Posicao no referencial PCPF
        h[i] = r[i] - simulacao.planeta.raio_equatorial
        # Forca propulsiva, massa e angulos
        ft[i], m[i], mu[i], epsl[i] = simulacao.foguete.modelo_propulsivo.propulsao_n_estagios(t[i], vetor_parametros,
                                                                                               simulacao.foguete.modelo_estrutural)
        # Parametros atmosfericos
        T[i], _, _, rho[i], _, Mach[i], _, _, Kn, _, _, R = simulacao.planeta.modelo_atmosferico.calcula(h[i], V[i],
                                                                                                         1.5, 10)

        simulacao.foguete.modelo_aerodinamico.atualizar_parametros(altitude=h[i], numero_de_knudsen=Kn,
                                                                   numero_de_mach=Mach[i],
                                                                   temperatura=T[i], constante_do_gas_ideal=R,
                                                                   velocidade=V[i])
        # Forcas aerodinamicas
        areas_de_referencia_para_calculo_do_arrasto, comprimento_caracteristico, fator_correcao = (
            simulacao.foguete.modelo_estrutural.calcula())
        D[i], _, _ = simulacao.foguete.modelo_aerodinamico.aerodinamica_multiplos_estagios(t[i], V[i],
                                                                                           areas_de_referencia_para_calculo_do_arrasto,
                                                                                           simulacao.foguete.modelo_propulsivo.tempos_de_separacao,
                                                                                           rho[i])
        # Pressao dinamica
        q[i] = 0.5 * rho[i] * V[i] ** 2
        # Coordenadas da velocidade inercial no referencial LVLH
        Vi[i], phii[i], Ai[i] = componentes_vel_relativa_para_inercial(V[i], phi[i], A[i], we, r[i], delta[i])
        # Longitude celeste
        longc[i] = converte_longitude_fixa_planeta_para_longitude_celeste(t[i], long[i], we, 0)
        # Energia especifica da orbita
        ee[i] = Vi[i] ** 2 / 2 - mut / r[i]
        # Posicao e velocidade inercial no referencial ICP
        rc0, vc0 = rvel_polar_para_rvel_retangular(Vi[i], Ai[i], phii[i], r[i], delta[i], longc[i])
        R0[i, :] = rc0.T
        rc01 = np.array([rc0[0], rc0[1], rc0[2]]).flatten()
        vc01 = np.array([vc0[0], vc0[1], vc0[2]]).flatten()
        a[i], e[i], in_[i], OM[i], om[i], _, tau[i] = calculos_orbitais.determina_parametros_orbitais(t[i], mut, rc01,
                                                                                                      vc01)

    # Analise de orbita
    # Altitude e velocidade inercial no fim da queima do terceiro estagio
    for i in range(N):
        if t[i] > simulacao.foguete.modelo_propulsivo.tempos_de_fim_de_queima[-1]:
            break

    ifq = i - 2  # Python's index is 0-based
    # Tempo do fim da queima do terceiro estagio
    tfq = t[ifq]
    # Velocidade inercial no fim da queima do terceiro estagio
    Vfq = Vi[ifq]
    # Altitude no fim da queima do terceiro estagio
    hfq = h[ifq]
    # Periodo da orbita obtida
    P = 2 * np.pi * np.sqrt((Requat + hfq) ** 3 / mut)
    print('*** Parametros da Orbita Obtida ***')
    print('Velocidade no momento da insercao orbital (km/s)', Vfq / 1e3)
    print('Altitude no momento da insercao orbital (km)', hfq)
    print('Distancia radial no momento da insercao orbital (km)', (hfq + Requat))
    print('Semi eixo maior (km)', a[ifq])
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
    print('Impulso de velocidade requerido para circularizacao da orbita (km/s)')
    DVgso = vgso - vagto
    print(DVgso / 1e3)
    print('Massa de propelente requerida para circularizacao da orbita (kg)')
    print('****** PARAMETROS DA ORBITA FINAL ******')
    print('Periodo (min)')
    P = 2 * np.pi * np.sqrt((a[-1] ** 3) / mut)
    print(P / 60)
    print('Semi eixo maior (km)', a[-1])
    print('Excentricidade', e[-1])
    print('Inclinacao (º)', in_[-1] * 180 / np.pi)

    plt.close('all')
    plt.figure(1)

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
    plt.ylabel('A (º)')

    plt.subplot(233)
    plt.plot(t, phi * 180 / np.pi, linewidth=2)
    plt.plot(tfq, phi[ifq - 1] * 180 / np.pi, '*')
    plt.grid(True)
    plt.axis('tight')
    plt.xlabel('t (s)')
    plt.ylabel('phi (º)')

    plt.subplot(234)
    plt.plot(t, h / 1e3, linewidth=2)
    plt.plot(t, np.full(N, hfq / 1e3), '--')
    plt.plot([tfq], [hfq / 1e3], '*')
    plt.grid(True)
    plt.axis('tight')
    plt.xlabel('t (s)')
    plt.ylabel('h (km)')
    plt.legend(['altitude', 'altitude no fim da queima do 3º estágio'])

    plt.subplot(235)
    plt.plot(t, delta * 180 / np.pi, linewidth=2)
    plt.grid(True)
    plt.axis('tight')
    plt.xlabel('t (s)')
    plt.ylabel('delta (º)')

    plt.figure(2)
    plt.subplot(221)
    plt.plot(t, Vi, linewidth=2)
    plt.plot(t, Vir, '--')
    plt.plot(t, np.full(N, Vfq), '-.')
    plt.grid(True)
    plt.xlabel('t (s)')
    plt.ylabel('V_i (m/s)')
    plt.legend(['Velocidade inercial', 'Velocidade de perigeu da órbita GTO requerida',
                'Velocidade no fim da queima do terceiro estágio'])

    plt.subplot(222)
    plt.plot(t, Ai * 180 / np.pi, linewidth=2)
    plt.grid(True)
    plt.axis('tight')
    plt.xlabel('t (s)')
    plt.ylabel('A_i (º)')

    plt.subplot(223)
    plt.plot(t, phii * 180 / np.pi, linewidth=2)
    plt.plot(tfq, phii[ifq - 1] * 180 / np.pi, '*')
    plt.grid(True)
    plt.axis('tight')
    plt.xlabel('t (s)')
    plt.ylabel('phi_i (º)')

    plt.subplot(224)
    plt.plot(t, longc * 180 / np.pi, linewidth=2)
    plt.grid(True)
    plt.axis('tight')
    plt.xlabel('t (s)')
    plt.ylabel('lambda (º)')

    plt.figure(3)
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
    plt.ylabel('\u03BC (º)')

    plt.subplot(224)
    plt.plot(t, epsl * 180 / np.pi, linewidth=2)
    plt.grid(True)
    plt.axis('tight')
    plt.xlabel('t (s)')
    plt.ylabel('\u03B5 (º)')

    plt.figure(4)
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
    plt.ylabel('T (ºC)')

    plt.subplot(326)
    plt.plot(t, rho, linewidth=2)
    plt.grid(True)
    plt.axis('tight')
    plt.xlabel('t (s)')
    plt.ylabel('rho (kg/m^3)')

    plt.figure(5)
    plt.subplot(311)
    plt.plot(t, ee, linewidth=2)
    plt.plot(t, eer, '--', linewidth=2)
    plt.plot(t, eegso, '--', linewidth=2)
    plt.grid(True)
    plt.xlabel('t (s)')
    plt.ylabel('\u03B5 (J/kg)')
    plt.legend(['Energia específica', 'Energia específica da órbita GTO requerida',
                'Energia específica da órbita GSO requerida'])

    plt.subplot(334)
    plt.plot(t, a / 1e3, linewidth=2)
    plt.plot(t, ar / 1e3, '--')
    plt.plot(t, Requat * np.ones([N, 1]) / 1e3, '-.')
    plt.grid(True)
    plt.xlabel('t (s)')
    plt.ylabel('a (km)')
    plt.legend(['Semi eixo maior', 'Semi eixo maior da órbita GTO requerida', 'Raio da Terra'])

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
    plt.ylabel('\u03C4 (s)')

    plt.subplot(337)
    plt.plot(t, OM * 180 / np.pi, linewidth=2)
    plt.grid(True)
    plt.axis('tight')
    plt.xlabel('t (s)')
    plt.ylabel('\u03A9 (º)')

    plt.subplot(338)
    plt.plot(t, in_ * 180 / np.pi, linewidth=2)
    plt.grid(True)
    plt.axis('tight')
    plt.xlabel('t (s)')
    plt.ylabel('i (º)')

    plt.subplot(339)
    plt.plot(t, om * 180 / np.pi, linewidth=2)
    plt.grid(True)
    plt.axis('tight')
    plt.xlabel('t (s)')
    plt.ylabel('\u03C9 (º)')

    plt.figure(7)
    ax = plt.axes(projection="3d")
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)
    r = Requat / 1e3

    x = r * np.outer(np.cos(u), np.sin(v))
    y = r * np.outer(np.sin(u), np.sin(v))
    z = r * np.outer(np.ones(np.size(u)), np.cos(v))

    ax.plot_surface(x, y, z, rstride=4, cstride=4)
    ax.plot3D(R0[:, 0] / 1e3, R0[:, 1] / 1e3, R0[:, 2] / 1e3, 'red')
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlabel('x (km)')
    ax.set_ylabel('y (km)')
    ax.set_zlabel('z (km)')

    #quero extrair a norma do primeiro vetor r0

    x,y,z = R0[0,:]

    r0 = np.array([x,y,z])
    plota_orbita(R0, simulacao.planeta.raio_equatorial, r0)
    plt.show()
