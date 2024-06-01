from Vrel2Vine  import Vrel2Vine
import numpy as np
from aquisicao_orbita_gso_na_mao import *
def parametros_manobra_adquire_gso(t, m, X):
    """
    Funcao para calcular os parametros da manobra mono impulsiva de aquisicao de orbita GSO
    Deve ser chamada no final da funcao de dinamica, pois precisa rodar ao
    longo do tempo da simulacao, verificando quando ocorre o apogeu da orbita
    GTO e determinando os parametros para imprimir o impulso de velocidade de
    circularizacao. Esta funcao nao entrega saidas, mas atualiza parametros
    por variaveis globais, que serao usados na funcao de propulsao
    Entradas:
    t (s): tempo
    m (kg): massa do foguete
    X: vetor de estado
    V = X[0] (m/s): modulo do vetor velocidade relativa com respeito ao
    planeta girante
    A = X[1] (rad): angulo de azimute do vetor velocidade relativa com
    respeito ao eixo z (que aponta para o norte) do sistema uen.
    phi = X[2] (rad): angulo de elevacao do vetor velocidade relativa com
    respeito ao horizonte local (plano yz do referencial uen)
    r = X[3] (m): distancia radial ate o centro do planeta
    delta = X[4] (rad): latitude com respeito ao plano equatorial do planeta
    lon = X[5] (rad): longitude planetaria
    """
    
    # Desmembra o vetor de estado
    V = X[0]
    A = X[1]
    phi = X[2]
    r = X[3]
    delta = X[4]
    # lon = X[5]

    # Vetor velocidade inercial
    Vi, phii, _ = Vrel2Vine(V, phi, A, we, r, delta)

    # Realizacao de uma sequencia de testes para verificar a ocorrencia do
    # apogeu da orbita GTO. Quando ele ocorre, determina os parametros da
    # manobra
    if r > 0.9 * agso:  # Soh inicia a verificacao quando esta perto do apogeu da orbita GTO    
        if not achouApogeu:  # Se nao achou o apogeu, entra na rotina de busca
            if np.sign(phii) != sinalPhii:  # Se o sinal for diferente, phii passou por zero, o foguete chegou no apogeu
                achouApogeu = True  # Encontrou o apogeu, ao setar essa variavel, soh vai entrar aqui uma vez
                ti[3] = t  # Tempo da segunda ignicao do motor do terceiro estagio
                
                # Calcula o tempo de duracao da queima para propiciar o DeltaV necessario
                DVgso = vgso - Vi  # Quando se nos passa nos testes acima, "vi" eh a velocidade de apogeu da GTO
                mp32 = (m * np.exp(DVgso / (Isp[2] * g)) - m) / np.exp(DVgso / (Isp[2] * g))  # Massa de propelente necessaria
                Tq32 = Tq3 * mp32 / mp3  # Duracao da queima necessaria
                
                # Verifica se o tempo eh maior que o maximo, se ocorrer, corrige
                if Tq31 + Tq32 > Tq3:
                    Tq32 = Tq3 - Tq31
                
                tq[3] = ti[3] + Tq32  # Tempo de fim de queima
                ts[2] = tq[3] + Ts3  # Tempo de separacao

    sinalPhii = np.sign(phii)  # Guarda o sinal de phi inercial para verificar mudanca na proxima iteracao
    
    return sinalPhii, achouApogeu, ti,tq,ts,Tq32