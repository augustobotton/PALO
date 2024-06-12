# Exemplo 5.7 do Tewari


# Importa o pacote feito em outra aula
# Importações de pacotes do Python
import numpy as np
# Inicia o exercício
# Constante gravitacional da Terra
mu=3.986e14; # m^3/s^2
# Geometria da órbita
a=6900 # km
a=a*1E3 # m
e=0.6
# Elementos orbitais da orbita inicial - dados
OMEGAi=120*np.pi/180;
omegai=25*np.pi/180;
ii=10*np.pi/180;
# Direcao do impulso de velocidade - dado
beta=100*np.pi/180;
alfa=(beta-np.pi/2)*2
# Distancia radial de apogeu
ra=(1+e)*a
# Velocidade de apogeu
va=np.sqrt(mu*(2/ra-1/a))
# Magnitude do impulso de velocidade
Dv=2*va*np.sin(alfa/2)
# Vetor velocidade no ponto de manobra no referencial perifocal inicial
vpi=np.array([0,-va,0])
# Vetor vetor posicao no ponto de manobra no referencial perifocal inicial
rpi=np.array([-ra,0,0])
# Vetor quantidade de movimento angular especifica da orbita inicial no
# referencial perifocal inicial
hpi=np.cross(rpi,vpi);
# Vetor impulso de velocidade responsavel pela manobra no referencial perifocal da
# primeira orbita
Dvp=Dv*(np.cos(beta)*vpi/va+np.sin(beta)*hpi/(ra*va))
# Vetor velocidade resultante da manobra (vetor velocidade de apogeu da
# orbita resultante) no referencial perifocal da órbita inicial
# orbita resultante) no referencial perifocal da órbita inicial
vpf=vpi+Dvp;
# Matriz de transformacao de coordenadas do referencial perifocal da órbita inicial
# para o celeste
C3O=np.array([[np.cos(OMEGAi), np.sin(OMEGAi), 0],[-np.sin(OMEGAi), np.cos(OMEGAi), 0],[0, 0, 1]])
C1i=np.array([[1, 0, 0],[0, np.cos(ii), np.sin(ii)],[0, -np.sin(ii), np.cos(ii)]])
C3o=np.array([[np.cos(omegai), np.sin(omegai), 0],[-np.sin(omegai), np.cos(omegai), 0],[0, 0, 1]])
Ccp=C3o@C1i@C3O
Cpc=np.transpose(Ccp)


# Transforma os vetores envolvidos para o referencial celeste (estão no referencial
# perifocal da órbita inicial)
vci=Cpc@vpi;rci=Cpc@rpi;hci=Cpc@hpi
vcf=Cpc@vpf
# Vetor posicao da orbita resultante no ponto de manobra escrito no
# referencial celeste. É o mesmo da órbita inicial
rcf=np.array(rci);
print(vcf);print(rcf)
## Elementos orbitais da orbita resultante
par_orb=det_orbita(0,rcf,vcf,mu)
print('Elementos orbitais da orbita resultante');
print('af = ',par_orb[0],'m');
print('ef = ',par_orb[1]);
print('OMEGAf = ',par_orb[3]*180/np.pi,' °');
print('if = ',par_orb[4]*180/np.pi,' °');
print('omegaf = ',par_orb[5]*180/np.pi,' °');
## Angulo entre as orbitas inicial e final
hcf=np.cross(rcf,vcf);
cosalpha=np.dot(hcf,hci)/(np.linalg.norm(hcf)*np.linalg.norm(hci))
alpha=np.arccos(cosalpha);
print('Angulo entre os planos orbitais inicial e final');
print('alfa = ',alpha*180/np.pi,' °');
