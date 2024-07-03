import numpy as np
from scipy.integrate import ode
import matplotlib.pyplot as plt

def example_10_02():
    # Constants
    global mu, R0, V0, t0, RE, J2
    mu = 398600  # Gravitational parameter (km^3/s^2)
    RE = 6378  # Earth's radius (km)
    J2 = 1082.63e-6  # J2 perturbation constant

    # Conversion factors
    hours = 3600  # Hours to seconds
    days = 24 * hours  # Days to seconds
    deg = np.pi / 180  # Degrees to radians

    # Initial orbital parameters (given)
    zp0 = 300  # Perigee altitude (km)
    za0 = 3062  # Apogee altitude (km)
    RA0 = 45 * deg  # Right ascension of the node (radians)
    i0 = 28 * deg  # Inclination (radians)
    w0 = 30 * deg  # Argument of perigee (radians)
    TA0 = 40 * deg  # True anomaly (radians)

    # Initial orbital parameters (inferred)
    rp0 = RE + zp0  # Perigee radius (km)
    ra0 = RE + za0  # Apogee radius (km)
    e0 = (ra0 - rp0) / (ra0 + rp0)  # Eccentricity
    a0 = (ra0 + rp0) / 2  # Semimajor axis (km)
    h0 = np.sqrt(rp0 * mu * (1 + e0))  # Angular momentum (km^2/s)
    T0 = 2 * np.pi / np.sqrt(mu) * a0**1.5  # Period (s)
    t0 = 0  # Initial time (s)
    tf = 2 * days  # Final time (s)

    # Store the initial orbital elements in the array coe0
    coe0 = [h0, e0, RA0, i0, w0, TA0]

    # Obtain the initial state vector from sv_from_coe function
    R0, V0 = sv_from_coe(coe0, mu)
    r0 = np.linalg.norm(R0)
    v0 = np.linalg.norm(V0)

    # Time step for Encke procedure
    del_t = T0 / 100

    # Begin the Encke integration
    t = t0
    tsave = [t0]
    y = np.hstack((R0, V0))
    del_y0 = np.zeros(6)
    t += del_t

    # Set up ODE solver
    solver = ode(rates).set_integrator('dopri5', max_step=del_t)
    solver.set_initial_value(del_y0, t0)

    # Integration loop
    while t <= tf + del_t / 2:
        solver.integrate(t)

        # Compute the osculating state vector at time t
        Rosc, Vosc = rv_from_r0v0(R0, V0, t - t0)

        # Rectify
        R0 = Rosc + solver.y[:3]
        V0 = Vosc + solver.y[3:]
        t0 = t

        # Prepare for next time step
        tsave.append(t)
        y = np.vstack((y, np.hstack((R0, V0))))
        t += del_t
        solver.set_initial_value(np.zeros(6), t)

    t = np.array(tsave)

    # Extract the orbital elements from the state vector
    n_times = len(t)
    h, e, RA, i, w, TA = [], [], [], [], [], []

    for j in range(n_times):
        R = y[j, :3]
        V = y[j, 3:]
        r = np.linalg.norm(R)
        v = np.linalg.norm(V)
        coe = coe_from_sv(R, V, mu)
        h.append(coe[0])
        e.append(coe[1])
        RA.append(coe[2])
        i.append(coe[3])
        w.append(coe[4])
        TA.append(coe[5])

    # Plot selected osculating elements
    plt.figure(1)
    plt.subplot(2, 1, 1)
    plt.plot(t / 3600, (np.array(RA) - RA0) / deg)
    plt.title('Variation of Right Ascension')
    plt.xlabel('hours')
    plt.ylabel(r'$\Delta \Omega$ (deg)')
    plt.grid(True, which='both')
    plt.tight_layout()

    plt.subplot(2, 1, 2)
    plt.plot(t / 3600, (np.array(w) - w0) / deg)
    plt.title('Variation of Argument of Perigee')
    plt.xlabel('hours')
    plt.ylabel(r'$\Delta \omega$ (deg)')
    plt.grid(True, which='both')
    plt.tight_layout()

    plt.figure(2)
    plt.subplot(3, 1, 1)
    plt.plot(t / 3600, np.array(h) - h0)
    plt.title('Variation of Angular Momentum')
    plt.xlabel('hours')
    plt.ylabel(r'$\Delta h$ (km$^2$/s)')
    plt.grid(True, which='both')
    plt.tight_layout()

    plt.subplot(3, 1, 2)
    plt.plot(t / 3600, np.array(e) - e0)
    plt.title('Variation of Eccentricity')
    plt.xlabel('hours')
    plt.ylabel(r'$\Delta e$')
    plt.grid(True, which='both')
    plt.tight_layout()

    plt.subplot(3, 1, 3)
    plt.plot(t / 3600, (np.array(i) - i0) / deg)
    plt.title('Variation of Inclination')
    plt.xlabel('hours')
    plt.ylabel(r'$\Delta i$ (deg)')
    plt.grid(True, which='both')
    plt.tight_layout()

    plt.show()

# Function to calculate the rates of Encke's deviation in position
def rates(t, f):
    global mu, R0, V0, t0, RE, J2

    # Unpack the state vector
    del_r = f[:3]
    del_v = f[3:]

    # Compute the state vector on the osculating orbit at time t (Equation 12.5)
    Rosc, Vosc = rv_from_r0v0(R0, V0, t - t0)

    # Calculate the components of the state vector on the perturbed orbit
    Rpp = Rosc + del_r
    Vpp = Vosc + del_v
    rosc = np.linalg.norm(Rosc)
    rpp = np.linalg.norm(Rpp)

    # Compute the J2 perturbing acceleration from Equation 12.30
    xx = Rpp[0]
    yy = Rpp[1]
    zz = Rpp[2]
    fac = 3 / 2 * J2 * (mu / rpp**2) * (RE / rpp)**2
    ap = -fac * np.array([(1 - 5 * (zz / rpp)**2) * (xx / rpp),
                          (1 - 5 * (zz / rpp)**2) * (yy / rpp),
                          (3 - 5 * (zz / rpp)**2) * (zz / rpp)])

    # Compute the total perturbing acceleration from Equation 12.7
    F = 1 - (rosc / rpp)**3
    del_a = -mu / rosc**3 * (del_r - F * Rpp) + ap

    # Time derivatives of the state vector
    dfdt = np.hstack((del_v, del_a))
    return dfdt

# Placeholder functions for sv_from_coe, coe_from_sv, and rv_from_r0v0
def sv_from_coe(coe, mu):
    # Placeholder implementation
    pass

def coe_from_sv(R, V, mu):
    # Placeholder implementation
    pass

def rv_from_r0v0(R0, V0, t):
    # Placeholder implementation
    pass

if __name__ == "__main__":
    example_10_02()
