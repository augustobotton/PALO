import numpy as np

Re = 6378.137*1000 #[km]
mu = 3.986e14 #[m^2/s^2]
h_i = 500*1000     #[km]
r_i = Re + h_i
a_f = 6900*1000
e_f = 0.6
r_pf = (1-e_f)*a_f
r_af = (1+e_f)*a_f
r_pt = r_i
r_at = r_af
a_t =(r_pt+r_at)/2
v_at =np.sqrt(mu*((2/r_at)-(1/a_t)))
v_pt =np.sqrt(mu*((2/r_pt)-(1/a_t)))
v_i =np.sqrt(mu*((2/r_i)-(1/r_i)))
v_f =np.sqrt(mu*((2/r_af)-(1/a_f)))

delta_v1 = v_pt - v_i
delta_v2 = v_f -  v_at
print(v_at)
print(v_pt)
print(v_i)
print(v_f)
print(delta_v1)
print(delta_v2)
print(abs(delta_v1)+abs(delta_v2))