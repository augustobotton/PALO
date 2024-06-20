# Gráficos
# Figure 1
import numpy as np
from matplotlib import pyplot as plt

plt.close('all')
plt.figure(1)

plt.subplot(231)
plt.plot(t, V, linewidth=2)
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('V (m/s)')

plt.subplot(232)
plt.plot(t, A * 180 / np.pi, linewidth=2)
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('A (º)')

plt.subplot(233)
plt.plot(t, phi * 180 / np.pi, linewidth=2)
plt.plot(tfq, phi[ifq - 1] * 180 / np.pi, '*')
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('phi (º)')

plt.subplot(234)
plt.plot(t, h / 1e3, linewidth=2)
plt.plot(t, hfq.T / 1e3, '--')
plt.plot(tfq, hfq[0][0] / 1e3, '*')
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('h (km)')
plt.legend(['altitude', 'altitude no fim da queima do 3º estágio'])

plt.subplot(235)
plt.plot(t, delta * 180 / np.pi, linewidth=2)
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('delta (º)')

plt.subplot(236)
plt.plot(t, lon * 180 / np.pi, linewidth=2)
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('l(º)')

# Figure 2
plt.figure(2)

plt.subplot(221)
plt.plot(t, Vi, linewidth=2)
plt.plot(t, Vir.T, '--')
plt.plot(t, Vfq.T, '-.')
plt.plot(tfq, Vfq[0][0], '*')
plt.grid(True)
plt.xlabel('t (s)')
plt.ylabel('V_i (m/s)')
plt.legend(['Velocidade inercial', 'Velocidade de perigeu da órbita GTO requerida',
            'Velocidade no fim da queima do terceiro estágio'])

plt.subplot(222)
plt.plot(t, Ai * 180 / np.pi, linewidth=2)
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('A_i (º)')

plt.subplot(223)
plt.plot(t, phii * 180 / np.pi, linewidth=2)
plt.plot(tfq, phii[ifq - 1] * 180 / np.pi, '*')
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('phi_i (º)')

plt.subplot(224)
plt.plot(t, longc * 180 / np.pi, linewidth=2)
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('lambda (º)')

# Figure 3
plt.figure(3)

plt.subplot(221)
plt.plot(t, ft, linewidth=2)
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('f_t (N)')

plt.subplot(222)
plt.plot(t, m, linewidth=2)
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('m (kg)')

plt.subplot(223)
plt.plot(t, mu * 180 / np.pi, linewidth=2)
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('\u03BC (º)')

plt.subplot(224)
plt.plot(t, epsl * 180 / np.pi, linewidth=2)
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('\u03B5 (º)')

# Figure 4
plt.figure(4)

plt.subplot(311)
plt.plot(t, D, linewidth=2)
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('D (N)')

plt.subplot(323)
plt.plot(t, q, linewidth=2)
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('q (N/m^2)')

plt.subplot(324)
plt.plot(t, M, linewidth=2)
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('M (-)')

plt.subplot(325)
plt.plot(t, T - 273.15, linewidth=2)
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('T (ºC)')

plt.subplot(326)
plt.plot(t, rho, linewidth=2)
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('rho (kg/m^3)')

# Figure 5
plt.figure(5)

plt.subplot(311)
plt.plot(t, ee, linewidth=2)
plt.plot(t, eer.T, '--', linewidth=2)
plt.plot(t, eegso.T, '--', linewidth=2)
plt.grid(True)
plt.xlabel('t (s)')
plt.ylabel('\u03B5 (J/kg)')
plt.legend(['Energia específica', 'Energia específica da órbita GTO requerida',
            'Energia específica da órbita GSO requerida'])

plt.subplot(334)
plt.plot(t, a / 1e3, linewidth=2)
plt.plot(t, ar.T / 1e3, '--')
plt.plot(t, raio_equatorial_terrestre * np.ones([N, 1]) / 1e3, '-.')
plt.grid(True)
plt.xlabel('t (s)')
plt.ylabel('a (km)')
plt.legend(['Semi eixo maior', 'Semi eixo maior da órbita GTO requerida', 'Raio da Terra'])

plt.subplot(335)
plt.plot(t, e, linewidth=2)
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('e (-)')

plt.subplot(336)
plt.plot(t, tau, linewidth=2)
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('\u03C4 (s)')

plt.subplot(337)
plt.plot(t, OM * 180 / np.pi, linewidth=2)
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('\u03A9 (º)')

plt.subplot(338)
plt.plot(t, inclinacao * 180 / np.pi, linewidth=2)
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('i (º)')

plt.subplot(339)
plt.plot(t, om * 180 / np.pi, linewidth=2)
plt.grid(True)
plt.axis('tight')
plt.xlabel('t (s)')
plt.ylabel('\u03C9 (º)')

# Figure 6
plt.figure(5)

traj = np.column_stack((delta, lon)) * 180 / np.pi
##desenha_mapa_trajetoria([delta0 * 180 / np.pi, lon0 * 180 / np.pi, h0], traj)
plt.show()

# Figure 7
fig7 = plt.figure(7)
ax = plt.axes(projection="3d")

u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
r = raio_equatorial_terrestre / 1e3

x = r * np.outer(np.cos(u), np.sin(v))
y = r * np.outer(np.sin(u), np.sin(v))
z = r * np.outer(np.ones(np.size(u)), np.cos(v))

ax.plot_surface(x, y, z, rstride=4, cstride=4)
ax.plot3D(R0[:, 0] / 1e3, R0[:, 1] / 1e3, R0[:, 2] / 1e3, 'red')
ax = plt.gca()
ax.set_aspect('equal', adjustable='box')
ax.set_xlabel('x (km)')
ax.set_ylabel('y (km)')
ax.set_zlabel('z (km)')

# Mostra os gráficos
plt.show()
