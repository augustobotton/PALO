
#Adaptado de André Luis da Silva
#
import numpy as np
#conda
def det_orbita(t0,rc0,vc0,mu):
    # Função para determinar parâmetros orbitais de órbita relativa de dois
    # corpos, a partir de uma observação de posição e velocidade no referencial
    # celeste no mesmo instante de tempo.
    # Entradas
    # t0: (s) tempo em que a observacao foi feita
    # rc0: (m ou km) vetor posição relativa de m2 com respeito a m1 escrito no
    # referencial celeste
    # vc0: (m/s ou km/s - unidades coerentes com o vetor posição) vetor
    # velocidade relativa escrito no referencial celeste. Deve ser tomado no
    # mesmo instante da medida da posição.
    # mu: (m**3/s**2 ou km**3/s**2 - unidades coerentes com a posição e velocidade)
    # parâmetro gravitacional padrão do corpo m1
    # Saídas
    # par_orb: vetor de parâmetros orbitais
    # a=par_orb(1): (m ou km - depende das unidades de entrada) semi eixo maior
    # da órbita
    # e=par_orb(2): (adimensional) excentricidade da órbita
    # tau=par_orb(3): (s) tempo de periastro
    # OMEGA=par_orb(4): (rad) longitude celeste do nodo ascendente
    # i=par_orb(5): (rad) inclinação
    # omega=par_orb(6): (rad) argumento de periastro
    ## Cálculos
    # Distância radial no instante observado
    r0=np.linalg.norm(rc0)
    # Vetor quantidade de movimento angular específica no referencial celeste
    hc=np.cross(rc0,vc0)
    # Vetor excentricidade no referencial celeste (eq. 2 da aula 9)
    ec=np.cross(vc0,hc)/mu-rc0/r0
    # Excentricidade da órbita
    e=np.linalg.norm(ec)
    # Módulo do vetor hc
    h=np.linalg.norm(hc)
    # Parâmetro da órbita
    p=h**2/mu
    # Semi eixo maior
    a=p/(1-e**2)
    # Vetor parâmetro no referencial celeste
    pc=p*np.cross(hc,ec)/(h*e)
    # Anomalia verdadeira no instante da observacao
    costheta=(p-r0)/(e*r0) # Eq. 5 da aula 8
    sintheta=np.dot(rc0,pc)/(r0*p) # Eq. 7 da aula 8
    theta0=np.arctan2(sintheta,costheta)
    # O tempo de periastro depende do tipo de órbita
    if (0<=e)and(e<1):
        tipo ='e' # Órbita elíptica
    elif e==1:
        tipo='p' # Órbita parabólica
    else:
        tipo='h' # Órbita hiperbólica

    # Tempo de periastro
    if tipo=='e': # Órbita elíptica
        # Movimento médio (Eq. 11 da aula 8)
        n=np.sqrt(mu/a**3)

        # Anomalia excêntrica (Eq. 10 da aula 8)
        E0=2*np.arctan(np.sqrt((1-e)/(1+e))*np.tan(theta0/2))
        # Tempo de periastro (Eq. 9 da aula 8)
        tau=t0-(E0-e*np.sin(E0))/n
    elif tipo=='p': # Órbita parabólica
        tau=-((np.tan(theta0/2)**3)+(3*np.tan(theta0/2)))/((mu/p**3)**(1/6)) # Eq. 12 da aula 9
    else: # Órbita hiperbólica
        # Anomalia hiperbólica (Eq. 14 da aula 8)
        n=np.sqrt(-mu/(a**3))
        H0=2*np.arctanh(np.sqrt((e-1)/(1+e))*np.tan(theta0/2))
        # Tempo de periastro (Eq. 13 da aula 8)
        tau=-(e*np.sinh(H0)-H0)/n


    # Linha dos nodos
    # Vetor unitário ao longo da linha dos nodos (no sistema celeste)
    ih=hc/h
    Kc=np.array([0,0,1])
    nc=np.cross(Kc,hc)/np.linalg.norm(np.cross(Kc,hc)) # (Eq. 15 da aula 8)
    # Longitude celeste do nodo ascendente (Eq. 18 da aula 8)
    OMEGA=np.arctan2(nc[1],nc[0])
    i = np.arccos(np.dot(ih,Kc))
    
    # Vetor unitário ao longo do vetor excentricidade (no referencial celeste)
    ie=ec/e
    # Argumento de periastro
    cosomega=np.dot(ie,nc) # (Eq. 22 da aula 8)
    sinomega=np.dot(ih,np.cross(nc,ie)) # (Eq. 23 da aula 8)
    omega=np.arctan2(sinomega,cosomega)
    ## Vetor de parâmetros de saída
    par_orb=np.empty(6)
    par_orb[0]=a;par_orb[1]=e;par_orb[2]=tau
    par_orb[3]=OMEGA;par_orb[4]=i;par_orb[5]=omega
    
    return par_orb