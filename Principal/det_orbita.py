import numpy as np

def det_orbita(t0, rc0, vc0, mu):
    # Função para determinar parâmetros orbitais a partir de uma observação de
    # posição e outra de velocidade, sendo as mesmas tomadas em relação ao
    # primário em um problema de dois corpos e escritas no referencial
    # celestial
    # Entradas
    # t0: tempo em que a observacao foi feita (s)
    # rc0: vetor posição em relação ao primário escrito no referencial
    # celeste (m ou km)
    # vc0: vetor velocidade em relação ao primário escrito no referencial
    # celeste (m/s ou km/s - unidades coerentes com o vetor posição). Deve
    # ser tomado no mesmo instante da medida da posição.
    # mu: parâmetro gravitacional padrão do corpo primário (m^3/s^2 ou km^3/s^2
    # - unidades coerentes com a posição e velocidade)
    # Saídas
    # par_orb: vetor de parâmetros orbitais
    # a=par_orb[0]: semi eixo maior da órbita (m ou km - depende das unidades
    # de entrada)
    # e=par_orb[1]: excentricidade da órbita (adimensional)
    # tau=par_orb[2]: tempo de perigeu da órbita (segundos)
    # OMEGA=par_orb[3]: ascenção reta do nodo ascendente (rad). Direção da
    # linha dos nodos da órbita dada em relação ao eixo X no plano XY do
    # referencial celeste
    # i=par_orb[4]: inclinação (rad). Inclinação da órbita dada em relação ao
    # plano XY do referencial celeste.
    # omega=par_orb[5]: argumento de perigeu (rad). Em relação à linha dos
    # nodos, medido no plano orbital.
    ## Cálculos
    # Distância radial ao primário no instante observado
    rc0 = rc0.T  # pra funcionar no np.cross
    vc0 = vc0.T

    r0 = np.linalg.norm(rc0);
    # Vetor quantidade de movimento angular específica no referencial celeste
    hc = np.cross(rc0, vc0);  # como rc0 e vc0 foram transpostos, hc vai ser 1x3 tambem
    # Vetor excentricidade no sistema celeste
    ec = np.cross(vc0, hc) / mu - rc0 / r0;  # voltando pra (3x1)
    # Excentricidade da órbita
    e = np.linalg.norm(ec);
    # Módulo do vetor hc
    h = np.linalg.norm(hc);

    # Parâmetro da órbita dada
    p = h ** 2 / mu;
    # Semi eixo maior

    a = p / (1 - e ** 2);

    # Vetor parâmetro no referencial celeste
    pc = p * np.cross(hc, ec) / (h * e);
    # Anomalia verdadeira
    costheta = (p - r0) / (e * r0);
    sintheta = np.dot(np.squeeze(rc0), np.squeeze(pc)) / (r0 * p);
    theta = np.arctan2(sintheta, costheta);
    # O tempo de perigeu depende do tipo de órbita
    if (0 <= e) and (e < 1):
        tipo = 'e';  # Órbita elíptica
    elif e == 1:
        tipo = 'p';  # Órbita parabólica
    else:
        tipo = 'h';  # Órbita hiperbólica

    # Tempo de perigeu
    if tipo == 'e':  # Órbita elíptica
        # Movimento médio
        n = np.sqrt(mu / a ** 3);
        # Anomalia excêntrica
        E = 2 * np.arctan(np.sqrt((1 - e) / (1 + e)) * np.tan(theta / 2));
        tau = t0 - (E - e * np.sin(E)) / n;
    elif tipo == 'p':  # Órbita parabólica
        tau = -((np.tan(theta / 2)) ** 3 + 3 * np.tan(theta / 2)) / (mu / p ** 3) ** (1 / 6);
    else:  # Órbita hiperbólica
        # Movimento médio hiperbólico
        n = np.sqrt(-mu / a ** 3);
        # Anomalia hiperbólica
        H = 2 * np.arctanh(np.sqrt((e - 1) / (1 + e)) * np.tan(theta / 2));
        tau = -(e * np.sinh(H) - H) / n;

    # Linha dos nodos
    # Vetor unitário ao longo do vetor h (no sistema celeste)
    ih = hc / h;
    # Vetor unitário ao longo da linha dos nodos (no sistema celeste)
    Kc = np.array([[0], [0], [1]]);
    nc = np.cross(Kc.T, ih) / np.linalg.norm(np.cross(Kc.T, ih));
    # Ascenção reta do nodo ascendente
    OMEGA = np.arctan2(nc[0][1], nc[0][0]);
    # Inclinação
    i = np.arccos(float(np.dot(ih, Kc)));
    # Vetor unitário ao longo do vetor excentricidade (no referencial
    # celeste)
    ie = ec / e;
    # Argumento de perigeu
    cosomega = np.dot(np.squeeze(ie), np.squeeze(nc));
    sinomega = np.dot(np.squeeze(ih), np.squeeze(np.cross(nc, ie)));
    omega = np.arctan2(sinomega, cosomega);
    ## Vetor de parâmetros de saída
    par_orb = [a, e, tau, OMEGA, i, omega];

    return par_orb