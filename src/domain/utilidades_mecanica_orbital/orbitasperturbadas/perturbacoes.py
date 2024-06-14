import numpy as np
# Função para cálculo de perturbações sobre órbitas Keplerianas
def perturb3Corpo(R,R31,mu3):
# Entradas:
# R: [m] - vetor 3x1. Posição da massa de prova m2 com respeito ao primário m1
# R31: [m] - vetor 3x1. Posição do terceiro corpo m3 com respeito ao primário m1
# Saída:
# ad: [m/s^2] - vetor 3x1. Aceleração perturbativa do terceiro corpo sobre m2
# Posição de m2 com respeito a m3
 R23=R-R31
 # Distâncias
 r23=np.linalg.norm(R23)
 r31=np.linalg.norm(R31)
 # Aceleração perturbativa
 ad=-mu3*(R23/r23**3+R31/r31**3)
 return ad