import numpy as np

def aceleracao_perturbativa_terceiro_corpo(
        posicao_m2: np.ndarray,
        posicao_m3: np.ndarray,
        parametro_gravitacional_m3: float
) -> np.ndarray:
    """
    Calcula a aceleração perturbativa de um terceiro corpo sobre uma massa de prova em órbita.

    :param posicao_m2: Vetor 3x1 da posição da massa de prova m2 com respeito ao primário m1 (m).
    :param posicao_m3: Vetor 3x1 da posição do terceiro corpo m3 com respeito ao primário m1 (m).
    :param parametro_gravitacional_m3: Parâmetro gravitacional do terceiro corpo m3 (m^3/s^2).
    :return: Vetor 3x1 da aceleração perturbativa do terceiro corpo sobre m2 (m/s^2).
    """
    # Posição de m2 com respeito a m3
    posicao_relativa_m2_m3 = posicao_m2 - posicao_m3

    # Distâncias
    distancia_m2_m3 = np.linalg.norm(posicao_relativa_m2_m3)
    distancia_m1_m3 = np.linalg.norm(posicao_m3)

    # Aceleração perturbativa
    aceleracao_perturbativa = -parametro_gravitacional_m3 * (
        posicao_relativa_m2_m3 / distancia_m2_m3 ** 3 + posicao_m3 / distancia_m1_m3 ** 3
    )

    return aceleracao_perturbativa
def dinPert3corpo(t,X, orbita):
    """
    Entradas:
    t: [s] - tempo
    X: vetor de estado, 6x1
    Saída
    xp: derivada de x
    """
    # Posição e velocidade
    R =np.array(X[0:3])
    V=np.array(X[3:6])
    #Distância
    r=np.linalg.norm(R)
    # Posição do terceiro corpo com respeito ao primário
    _,R31,_=propagacao.propagaEliptica(t,mu,orbita)
    # Aceleração perturbativa
    ad= aceleracao_perturbativa_terceiro_corpo(R, R31, mu3)
    # Cinemática
    Rp=V
    # Dinâmica
    Vp=-mu*R/r**3+ad
    #Saída
    Xp=np.concatenate((Rp,Vp))
    return Xp
