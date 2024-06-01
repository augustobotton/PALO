from principal_GSO import *

def parametros_manobra_adquire_gso(t, m, X):
    # Desmembra o vetor de estado
    global we, agso, sinalPhii, achouApogeu, ti, tq, ts, Tq3, Tq31, Tq32, Ts3, vgso, Isp, g, mp3
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