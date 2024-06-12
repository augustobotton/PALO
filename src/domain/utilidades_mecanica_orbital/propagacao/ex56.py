import numpy as np

Re = 6378.137*1000 #[km]
mu = 3.986e14 #[m^2/s^2]
h_i = 500*1000     #[km]
i_i = 10      #[grau]
h_pf = 200*1000    #[km]
h_af = 700*1000    #[km]
i_f = 5       #[grau]

r_i  = Re + h_i
r_pf = Re + h_pf
r_af = Re + h_af

a_f = (r_pf+r_af)/2
v_f = np.sqrt(mu*((2/r_i)-(1/a_f)))
print("Velocidade da Nova Órbita:", v_f)
e_f = r_af/a_f -1
print("Excentricidade Final" ,e_f)
p_f = r_af*(1-e_f)
print("Parâmetro P",p_f)
h_f = np.sqrt(p_f*mu)
print("Altitude Final", h_f)

cosphi =h_f/(r_i*v_f)
alpha=np.arccos(cosphi)
print("Ângulo de Trajetória", np.rad2deg(alpha) )
v_i = np.sqrt(mu*((2/r_i)-(1/r_i)))
print("Velocidade Inicial Órbita Circular:", v_i)
delta_v = np.sqrt(v_i**2+v_f**2-(2*v_i*v_f*cosphi))
print("Primeira Variação da Magnitude da Velocidade: ",delta_v)

betha = np.arctan((v_f*np.sin(alpha))/(v_f*np.cos(alpha)-v_i))
print("Beta 1:",np.rad2deg(betha))


delta_inclinacao = i_i-i_f
print("Variação Requerida na inclinação da Órbita:", delta_inclinacao)
vi_aposmanobra = np.sqrt(mu*((2/r_af)-(1/a_f)))
print("Velocidade Apogeu após manobra:", vi_aposmanobra)

deltav_segundo = 2*vi_aposmanobra*np.sin(np.deg2rad(delta_inclinacao)/2)
print("Segundo impulso de velocidade:",deltav_segundo)
betha2 = np.rad2deg(np.pi/2)+delta_inclinacao/2
print("Beta 2",betha2)