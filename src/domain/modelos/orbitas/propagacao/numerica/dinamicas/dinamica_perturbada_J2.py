from src.domain.modelos.orbitas.utilidades.rv_from_r0v0 import rv_from_r0v0
import numpy as np

def dinamica_perturbada_J2(t, f, *args):
    mu = args[0]
    R0 = args[1]
    V0 = args[2]
    t0 = args[3]
    RE = args[4]
    J2 = args[5]

    # Unpack the state vector
    del_r = f[:3]
    del_v = f[3:6]
    # Compute the state vector on the osculating orbit at time t
    Rosc, Vosc = rv_from_r0v0(R0, V0, t - t0, mu)
    # Calculate the components of the state vector on the perturbed orbit
    Rpp = Rosc + del_r
    Vpp = Vosc + del_v
    rosc = np.linalg.norm(Rosc)
    rpp = np.linalg.norm(Rpp)

    # Compute the J2 perturbing acceleration
    xx = Rpp[0]
    yy = Rpp[1]
    zz = Rpp[2]
    fac = 1.5 * J2 * (mu / rpp ** 2) * (RE / rpp) ** 2
    ap = -fac * np.array([
        (1 - 5 * (zz / rpp) ** 2) * (xx / rpp),
        (1 - 5 * (zz / rpp) ** 2) * (yy / rpp),
        (3 - 5 * (zz / rpp) ** 2) * (zz / rpp)
    ])

    # Compute the total perturbing acceleration
    F = 1 - (rosc / rpp) ** 3
    del_a = (-mu / rosc ** 3) * (del_r - F * Rpp) + ap

    return np.concatenate([del_v, del_a])