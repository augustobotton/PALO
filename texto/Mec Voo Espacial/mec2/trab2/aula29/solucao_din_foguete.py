import numpy as np
from aerodinamica_N_estagios import aerodinamica_N_estagios
from atm_padrao import atm_padrao
from propulsao_N_estagios import propulsao_N_estagios
from grav_axisimetrico import grav_axissimetrico
from parametros_manobra_adquire_gso import parametros_manobra_adquire_gso
import aquisicao_orbita_gso_na_mao as ent

def solucao_din_foguete(t, X):
        # Vetor de estado
                
        V, A, phi, r, delta = X[0], X[1], X[2], X[3], X[4]
       
        if V < 0:
           V = 0  # Evita velocidade negativa
            
        # Funcao para calculo da massa e tracao em funcao do tempo
        ft, m, mu, epsl = propulsao_N_estagios(t, V, A, phi, r, delta)
        
        # Funcao para calculo do modelo atmosferico
        h = r - ent.Requat
        
        T, Tm, p, rho, ainf, Mach, _, Pr, Kn, d, Reynolds,R = atm_padrao(h, V, ent.lc, ent.dT)
        #print("mach:", Mach,"Tempo:",t,"altitude: ",h,"velocidade:",V)
        # Calculo do modelo aerodinamico
        D, fy, L = aerodinamica_N_estagios(t, V, h, Mach, Kn, T, rho, R,ent.fc,ent.ts,ent.Sr)
        
        # Calculo da gravidade
        gc, gd = grav_axissimetrico(r, delta)
        # Equacoes de cinematica de translacao
        rp = V * np.sin(phi)
        deltap = (V / r) * np.cos(phi) * np.cos(A)
        lonp = (V * np.cos(phi) * np.sin(A)) / (r * np.cos(delta))
        # Equacoes de dinamica de translacao

        Vp = (1 / m) * (ft * np.cos(epsl) * np.cos(mu) - D - m * gc * np.sin(phi) + m * gd * np.cos(phi) * np.cos(A) - m * ent.we**2 * r * np.cos(delta) * (np.cos(phi) * np.cos(A) * np.sin(delta) - np.sin(phi) * np.cos(delta)))
        Ap = (1 / (m * V * np.cos(phi))) * (m * (V**2 / r) * np.cos(phi)**2 * np.sin(A) * np.tan(delta) + ft * np.sin(mu) + fy - m * gd * np.sin(A) + m * ent.we**2 * r * np.sin(A) * np.sin(delta) * np.cos(delta) - 2 * m * ent.we * V * (np.sin(phi) * np.cos(A) * np.cos(delta) - np.cos(phi) * np.sin(delta)))
        phip = (1 / (m * V)) * (m * (V**2 / r) * np.cos(phi) + ft * np.sin(epsl) * np.cos(mu) + L - m * gc * np.cos(phi) - m * gd * np.sin(phi) * np.cos(A) + m * ent.we**2 * r * np.cos(delta) * (np.sin(phi) * np.cos(A) * np.sin(delta) + np.cos(phi) * np.cos(delta)) + 2 * m * ent.we * V * np.sin(A) * np.cos(delta))
        # Saturacao da altitude
        if h < 0:  # Altitude negativa nao e permitida
            # Mantem as derivadas nulas
            rp = 0
            deltap = 0
            lonp = 0
            Vp = 0
            Ap = 0
            phip = 0
        # Modela o trilho de lancamento
        H = h - ent.h0  # Altura
        if (H <= ent.l_trilho) and t <= 10:  # Verifica se a altura eh menor que l_trilho nos primeiros segundos da simulacao
            Ap = 0
            phip = 0  # Anula as derivadas dos angulos de orientacao da velocidade  
        # Derivada do vetor de estado
        
        parametros_manobra_adquire_gso(t, m, X,)

        Xp = [Vp, Ap, phip, rp, deltap, lonp]
    
        return Xp
