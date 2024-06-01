import numpy as np

def det_orbita(t0, rc0, vc0, mu):
    # Calculations
    # Radial distance to the primary at the observed instant
    r0 = np.linalg.norm(rc0)
    # Vector of specific angular momentum in the celestial reference frame
    hc = np.cross(rc0, vc0)
    # Eccentricity vector in the celestial system
    ec = np.cross(vc0, hc) / mu - rc0 / r0
    # Orbit eccentricity
    e = np.linalg.norm(ec)
    # Module of the vector hc
    h = np.linalg.norm(hc)
    # Parameter of the given orbit
    p = h ** 2 / mu
    # Major semi-axis
    a = p / (1 - e ** 2)
    # Parameter vector in the celestial reference frame
    pc = p * np.cross(hc, ec) / (h * e)
    # True anomaly
    costheta = (p - r0) / (e * r0)
    sintheta = np.dot(rc0, pc) / (r0 * p)
    theta = np.arctan2(sintheta, costheta)
    # Perigee time depends on the type of orbit
    if (0 <= e) and (e < 1):
        tipo = 'e'   # Elliptical orbit
    elif e == 1:
        tipo = 'p'   # Parabolic orbit
    else:
        tipo = 'h'   # Hyperbolic orbit

    # Perigee time
    if tipo == 'e':   # Elliptical orbit
        # Mean motion
        n = np.sqrt(mu / a ** 3)
        # Eccentric anomaly
        E = 2 * np.arctan(np.sqrt((1 - e) / (1 + e)) * np.tan(theta / 2))
        tau = t0 - (E - e * np.sin(E)) / n
    elif tipo == 'p':   # Parabolic orbit
        tau = -((np.tan(theta / 2)) ** 3 + 3 * np.tan(theta / 2)) / (mu / p ** 3) ** (1 / 6)
    else:   # Hyperbolic orbit
        # Hyperbolic mean motion
        n = np.sqrt(-mu / a ** 3)
        # Hyperbolic anomaly
        H = 2 * np.arctanh(np.sqrt((e - 1) / (1 + e)) * np.tan(theta / 2))
        tau = -(e * np.sinh(H) - H) / n

    # Line of nodes
    # Unit vector along the h vector (in the celestial system)
    ih = hc / h
    # Unit vector along the line of the nodes (in the celestial system)
    Kc = [0, 0, 1]
    nc = np.cross(Kc, ih) / np.linalg.norm(np.cross(Kc, ih))
    # Right ascension of the ascending node
    OMEGA = np.arctan2(nc[1], nc[0])
    # Inclination
    i = np.arccos(np.dot(ih, Kc))
    # Unit vector along the eccentricity vector (in the celestial reference frame)
    ie = ec / e
    # Argument of perigee
    cosomega = np.dot(ie, nc)
    sinomega = np.dot(ih, np.cross(nc, ie))
    omega = np.arctan2(sinomega, cosomega)

    # Output parameters vector
    par_orb = [a, e, tau, OMEGA, i, omega]
    return par_orb