import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks

RE=6378

def calcular_magnitudes(y):
    """
    Calcula as magnitudes dos vetores posição e velocidade.

    Parâmetros:
    y (ndarray): Matriz contendo vetores posição e velocidade.

    Retorna:
    tuple: Magnitudes dos vetores posição e velocidade.
    """
    r = np.linalg.norm(y[:, :3], axis=1)
    v = np.linalg.norm(y[:, 3:6], axis=1)
    return r, v

def encontrar_extremos(r, t, RE):
    """
    Encontra os valores extremos de altitude e suas correspondentes velocidades e tempos.

    Parâmetros:
    r (ndarray): Array vetor Posição normalizado.
    t (ndarray): Vetor de tempos.
    RE (float): Raio da Terra.

    Retorna:
    tuple: Altitudes e tempos dos máximos e mínimos locais, além de outras informações relevantes.
    """
    rmax = np.max(r)
    rmin = np.min(r)
    imax = np.argmax(r)
    imin = np.argmin(r)
    altitude = r - RE

    imax_local, _ = find_peaks(altitude)
    imin_local, _ = find_peaks(-altitude)
    maxima = np.column_stack((t[imax_local], altitude[imax_local]))
    minima = np.column_stack((t[imin_local], altitude[imin_local]))

    apogee = maxima[maxima[:, 0].argsort()]
    perigee = minima[minima[:, 0].argsort()]

    return rmax, rmin, imax, imin, apogee, perigee

def plotar_apogeu_perigeu(apogee, perigee):
    """
    Plota as altitudes máximas e mínimas ao longo do tempo.

    Parâmetros:
    apogee (ndarray): Tempos e altitudes dos apogeus.
    perigee (ndarray): Tempos e altitudes dos perigeus.
    """
    plt.figure(1)

    # Convertendo o tempo de segundos para dias
    tempo_apogee_horas = apogee[:, 0] / (3600*24)
    tempo_perigee_horas = perigee[:, 0] / (3600*24)

    plt.plot(tempo_apogee_horas, apogee[:, 1], 'b', linewidth=2, label='Apogeu')
    plt.plot(tempo_perigee_horas, perigee[:, 1], 'r', linewidth=2, label='Perigeu')
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.xlabel('Tempo (Dias)')
    plt.ylabel('Altitude (km)')
    plt.ylim([0, 1000])
    plt.title('Altitude na Trajetória Orbital', fontsize=14, fontweight='bold')
    plt.legend()
    plt.tight_layout()
    plt.show()

def imprimir_resultados(r0, v0, t, r, rmax, rmin, imax, imin, apogee, perigee, y):
    """
    Imprime os resultados da simulação da órbita terrestre.

    Parâmetros:
    r0 (ndarray): Vetor posição inicial.
    v0 (ndarray): Vetor velocidade inicial.
    t (ndarray): Vetor de tempos.
    r (ndarray): Magnitudes dos vetores posição.
    rmax (float): Altitude máxima.
    rmin (float): Altitude mínima.
    imax (int): Índice da altitude máxima.
    imin (int): Índice da altitude mínima.
    apogee (ndarray): Tempos e altitudes dos apogeus.
    perigee (ndarray): Tempos e altitudes dos perigeus.
    y (ndarray): Matriz contendo vetores posição e velocidade.
    """
    print("\nÓrbita Terrestre")
    print("A posição inicial é", r0)
    print("Magnitude =", np.linalg.norm(r0), "km")
    print("A velocidade inicial é", v0)
    print("Magnitude =", np.linalg.norm(v0), "km/s")
    print("Tempo inicial = {:.2f} h. Tempo final = {:.2f} h.".format(t[0] / 3600, t[-1] / 3600))
    print("A altitude mínima é {:.2f} km no tempo = {:.2f} h.".format(rmin - RE, t[imin] / 3600))
    print("A velocidade nesse ponto é {:.2f} km/s.".format(np.linalg.norm(y[imin, 3:6])))
    print("A altitude máxima é {:.2f} km no tempo = {:.2f} h.".format(rmax - RE, t[imax] / 3600))
    print("A velocidade nesse ponto é {:.2f} km/s.".format(np.linalg.norm(y[imax, 3:6])))
    r_final = y[-1, :3]
    v_final = y[-1, 3:6]
    print("\nVetor posição no tempo final:", r_final)
    print("Vetor velocidade no tempo final:", v_final)

def obter_vetor_estado_apoastro(t, y, apogee):
    """
    Obtém o vetor de estado (posição e velocidade) no apoastro.

    Parâmetros:
    t (ndarray): Vetor de tempos.
    y (ndarray): Matriz contendo vetores posição e velocidade.
    apogee (ndarray): Tempos e altitudes dos apogeus.

    Retorna:
    tuple: Vetor posição e vetor velocidade no apoastro.
    """
    if len(apogee) == 0:
        return None, None  # Nenhum apoastro encontrado

    t_apoastro = apogee[0, 0]  # Tempo do primeiro apoastro
    idx_apoastro = np.where(t == t_apoastro)[0][0]
    pos_apoastro = y[idx_apoastro, :3]
    vel_apoastro = y[idx_apoastro, 3:6]

    return pos_apoastro, vel_apoastro

def constroi_resultados(t, y, RE, r0, v0):
    """
    Constrói e exibe os resultados da simulação da órbita terrestre.

    Parâmetros:
    t (ndarray): Vetor de tempos.
    y (ndarray): Matriz contendo vetores posição e velocidade.
    RE (float): Raio da Terra.
    r0 (ndarray): Vetor posição inicial.
    v0 (ndarray): Vetor velocidade inicial.
    """
    y = np.asarray(y)
    if y.ndim == 1:
        raise ValueError("A matriz 'y' deve ser 2-dimensional.")

    r, v = calcular_magnitudes(y)
    rmax, rmin, imax, imin, apogee, perigee = encontrar_extremos(r, t, RE)
    plotar_apogeu_perigeu(apogee, perigee)
    imprimir_resultados(r0, v0, t, r, rmax, rmin, imax, imin, apogee, perigee, y)

    pos_apoastro, vel_apoastro = obter_vetor_estado_apoastro(t, y, apogee)
    if pos_apoastro is not None and vel_apoastro is not None:
        print("\nVetor posição no apoastro:", pos_apoastro)
        print("Vetor velocidade no apoastro:", vel_apoastro)





