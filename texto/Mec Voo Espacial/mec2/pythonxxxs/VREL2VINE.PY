import math as m

def vRel2vIne(VR, PHIR, AR, WE, R, DT):
    #   ENTRADAS
    #
    #   VARIÁVEL    UNIDADE DESCRIÇÃO
    #   VR          m/s     Velocidade relativa
    #   PHIR        rad     Inclinação do vetor velocidade relativa
    #   AR          rad     Azimute do vetor velocidade relativa
    #   WE          rad/s   Velocidade de rotação do referencial girante
    #   R           m       Distância radial até a origem do ref. inercial
    #   dt
    #
    #   SAÍDAS
    #
    #   VARIÁVEL    UNIDADE DESCRIÇÃO
    #   VI          m/s     Magnitude da velocidade ref. inercial
    #   PHII        rad     Ângulo de elevação do vetor velocidade inercial
    #   AI          rad     Ângulo de azimute do vetor velocidade inercial
    
    AI = m.atan2(   VR*m.cos(PHIR)*m.sin(PHIR)  \
                    + WE*R*m.cos(DT),           \
                                                \
                    VR*m.cos(PHIR)*m.cos(AR)    )
    
    if(AI < 0):
        AI = AI + 2*m.pi

    VI = m.sqrt(    VR**2 + 2*VR*m.cos(PHIR)*m.sin(AR)*R*WE*m.cos(DT)   \
                    + (R**2)*(WE**2)*(m.cos(DT)**2)                     )
    
    
    PHII = m.atan2(m.sin(PHIR)*m.cos(AI), m.cos(PHIR)*m.cos(AR))

    if(abs(PHII) > m.pi/2.0):
            
            if(PHII < m.pi/2.0):
                    PHII = PHII + m.pi
            
            elif(PHII > m.pi/2.0):
                    PHII = PHII - m.pi

    return VI, PHII, AI
