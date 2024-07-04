import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Definindo variáveis globais
G = 6.6742e-20  # Universal gravitational constant (km^3/kg/s^2)
m1 = 5.974e24  # Planet mass (kg)
m2 = 1000  # Spacecraft mass (kg)
R = 6378  # Planet radius (km)
mu = G * (m1 + m2)





def orbit(ti, tf, r0, v0):
    hours = 3600
    r0 = [8000, 0, 6000]
    v0 = [0, 7, 0]
    ti = 0
    tf = 4 * hours

    y0 = np.array(r0 + v0)
    opts = {'rtol': 1e-8, 'atol': 1e-8}
    sol = solve_ivp(dinamica_lambert, [ti, tf], y0, **opts)

    t = sol.t
    y = sol.y.T

    r = np.linalg.norm(y[:, :3], axis=1)
    rmax = np.max(r)
    rmin = np.min(r)
    imax = np.argmax(r)
    imin = np.argmin(r)

    v_at_rmax = np.linalg.norm(y[imax, 3:6])
    v_at_rmin = np.linalg.norm(y[imin, 3:6])

    print("\nEarth Orbit")
    print("The initial position is", r0)
    print("Magnitude =", np.linalg.norm(r0), "km")
    print("The initial velocity is", v0)
    print("Magnitude =", np.linalg.norm(v0), "km/s")
    print("Initial time = {:.2f} h. Final time = {:.2f} h.".format(ti / hours, tf / hours))
    print("The minimum altitude is {:.2f} km at time = {:.2f} h.".format(rmin - R, t[imin] / hours))
    print("The speed at that point is {:.2f} km/s.".format(v_at_rmin))
    print("The maximum altitude is {:.2f} km at time = {:.2f} h.".format(rmax - R, t[imax] / hours))
    print("The speed at that point is {:.2f} km/s.".format(v_at_rmax))

    plot_orbit(y, R, r0)


def plot_orbit(y, R, r0):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    u, v = np.mgrid[0:2 * np.pi:100j, 0:np.pi:50j]
    x = R * np.cos(u) * np.sin(v)
    y_sphere = R * np.sin(u) * np.sin(v)
    z = R * np.cos(v)
    ax.plot_surface(x, y_sphere, z, color='grey', alpha=0.5)  # Desenha o planeta

    ax.plot(y[:, 0], y[:, 1], y[:, 2], 'k')  # Plota a trajetória da órbita
    ax.plot([0, r0[0]], [0, r0[1]], [0, r0[2]], 'r')  # Linha radial do ponto inicial
    ax.text(y[0, 0], y[0, 1], y[0, 2], 'o', color='red')  # Marca o ponto inicial
    ax.text(y[-1, 0], y[-1, 1], y[-1, 2], 'f', color='blue')  # Marca o ponto final

    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')
    ax.set_zlabel('Z (km)')
    ax.set_aspect('auto')
    plt.show()


if __name__ == "__main__":
    orbit()
