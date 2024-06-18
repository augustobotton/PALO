# Exemplo 5.7 do Tewari

# Importações de pacotes do Python
import numpy as np

# Inicia o exercício
# Constante gravitacional da Terra
mu=3.986e14; # m^3/s^2

# Geometria da órbita
a=6900*1E3 # m
e=0.6
# Elementos orbitais da orbita inicial - dados
OMEGAi=120*np.pi/180;
omegai=25*np.pi/180;
ii=10*np.pi/180;



# Direcao do impulso de velocidade - dado
beta=100*np.pi/180;
alfa=(beta-np.pi/2)*2

# Distancia radial de apogeu # adicinar funcao em orbita
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

Cpc = matriz_rotacao_orbital_inercial(OMEGAi, omegai, ii)
# Transforma os vetores envolvidos para o referencial celeste (estão no referencial
# perifocal da órbita inicial)
vci=Cpc@vpi;
rci=Cpc@rpi;
hci=Cpc@hpi
vcf=Cpc@vpf
# Vetor posicao da orbita resultante no ponto de manobra escrito no
# referencial celeste. É o mesmo da órbita inicial
rcf=np.array(rci);
print(vcf);print(rcf)
## Elementos orbitais da orbita resultante

par_orb=determina_parametros_orbitais(0,mu,rcf,vcf)
print('Elementos orbitais da orbita resultante');
print('af = ',par_orb.semi_eixo_maior,'m');
print('ef = ',par_orb.excentricidade);
print('OMEGAf = ',np.rad2deg(par_orb.raan),' °');
print('if = ',np.rad2deg(par_orb.inclinacao),' °');
print('omegaf = ',np.rad2deg(par_orb.arg_periastro),' °');
## Angulo entre as orbitas inicial e final
hcf=np.cross(rcf,vcf);
cosalpha=np.dot(hcf,hci)/(np.linalg.norm(hcf)*np.linalg.norm(hci))
alpha=np.arccos(cosalpha);
print('Angulo entre os planos orbitais inicial e final');
print('alfa = ',alpha*180/np.pi,' °');



