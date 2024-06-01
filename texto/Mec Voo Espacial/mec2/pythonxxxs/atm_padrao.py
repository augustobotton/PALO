import numpy as np

def atm_padrao(h, v, lc, dT):

    # Constants
    hi = np.array([0, 11.0191, 20.0631, 32.1619, 47.3501, 51.4125, 71.8020, 86, 100, 110, 120, 150,
                   160, 170, 190, 230, 300, 400, 500, 600, 700]) * 1e3
    Ti = np.array([288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65, 186.946, 210.65, 260.65, 
                   360.65, 960.65, 1110.65, 1210.65, 1350.65, 1550.65, 1830.65, 2160.65, 2420.65, 
                   2590.65, 2700.65])
    Ri = np.array([287.0, 287.0, 287.0, 287.0, 287.0, 287.0, 287.02, 287.02, 287.84, 291.06, 308.79, 
                   311.80, 313.69, 321.57, 336.68, 366.84, 416.88, 463.36, 493.63, 514.08, 514.08])
    a = np.array([-6.5, 0.0, 1.0, 2.8, 0.0, -2.8, -2.0, 1.693, 5.0, 10.0, 20.0, 15.0, 10.0, 7.0, 5.0, 
                  4.0, 3.3, 2.6, 1.7, 1.1, 0.0]) * 1e-3
    Mi = np.array([28.9644, 28.9644, 28.9644, 28.9644, 28.9644, 28.9644, 28.9644, 28.9644, 
                   28.88, 28.56, 28.07, 26.92, 26.66, 26.4, 25.85, 24.7, 22.66, 19.94, 17.94, 16.84, 
                   16.17])

    M0 = Mi[0]  # Molecular weight at sea level
    g0 = 9.80665  # Gravitational acceleration at sea level (m/s^2)
    Na = 6.0220978e23  # Avogadro's number
    sigma = 3.65e-10  # Collision diameter for air (m)
    m0 = 28.964e-3  # Molar mass of air at sea level (kg/Mol)
    P0 = 1.01325e5  # Standard pressure at sea level (N/m^2)
    Re = 6378.14e3  # Earth's mean radius (m)
    gamma = 1.405  # Specific heat ratio at sea level

    # Beta number associated with the mean radial distance from sea level
    beta = 2 / Re

    # Identify the layer to which the altitude belongs
    i = 1
    if h < 0:
        i = 1
        h = 0
    elif h > 2000e3:
        i = 21
        h = 2000e3
    else:
        for i in range(21):
            if (i == 20) or ((h >= hi[i]) and (h < hi[i+1])):
                break

    # Pressure at the beginning of layer i
    Pi = P0  # Initial pressure of the layer - initialized with the sea level value
    px = P0  # Auxiliary variable to store the initial pressure value in the previous layer
    for j in range(1, i):
        if a[j-1] != 0:  # Check if the layer is not isothermal
            A = 1 + (a[j-1] * (hi[j] - hi[j-1])) / Ti[j-1]
            B = -(g0 / (Ri[j] * a[j-1])) * (1 + beta * (Ti[j-1] / a[j-1] - hi[j-1]))
            C = (g0 * beta / (Ri[j] * a[j-1])) * (hi[j] - hi[j-1])
            Pi = px * (A ** B) * np.exp(C)
            px = Pi
        else:
            Pi = px * np.exp(-(g0 / (Ri[j] * Ti[j-1])) * (hi[j] - hi[j-1]) * 
                             (1 - (beta / 2) * (hi[j] - hi[j-1])))
            px = Pi

    Tm = Ti[i] + a[i] * (h - hi[i])
    # Correct for the chosen delta T value
    Tm = Tm + dT
    # Interpolate the value of R and M
    if i >= 20:
        R = Ri[i]
        M = Mi[i]
    else:
        R = Ri[i] + ((Ri[i+1] - Ri[i]) / (hi[i+1] - hi[i])) * (h - hi[i])
        M = Mi[i] + ((Mi[i+1] - Mi[i]) / (hi[i+1] - hi[i])) * (h - hi[i])
        
    # Standard (kinetic) temperature - function output
    T = (M / M0) * Tm
    # Pressure
    if a[i] != 0:  # Check if the layer is not isothermal
        # Calculate pressure
        A = 1 + (a[i] * (h - hi[i])) / Ti[i]
        B = -(g0 / (R * a[i])) * (1 + beta * (Ti[i] / a[i] - hi[i]))
        C = (g0 * beta / (R * a[i])) * (h - hi[i])
        p = Pi * (A ** B) * np.exp(C)
    else:
        # Calculate pressure using the isothermal model
        p = Pi * np.exp(-(g0 / (R * Ti[i])) * (h - hi[i]) * (1 - (beta / 2) * (h - hi[i])))
        
    # Calculate density
    rho = p / (R * Tm)
    # Speed of sound
    ainf = np.sqrt(gamma * R * Tm)
    # Mach number
    M = v / ainf
    # Dynamic viscosity coefficient
    mu = 1.458e-6 * (Tm ** 1.5) / (Tm + 110.4)
    # Prandtl number
    cp = R * gamma / (gamma - 1)
    kT = (2.64638e-3 * (Tm ** 1.5)) / (Tm + 245.4 * (10 ** (-12 / Tm)))
    Pr = mu * cp / kT
    # Knudsen number
    lam = m0 / (np.sqrt(2) * np.pi * sigma ** 2 * rho * Na)
    Kn = lam / lc
    # Flow regime parameter
    if Kn >= 10:
        d = 1
    elif Kn <= 0.01:
        d = 2
    else:
        d = 3
        
    # Reynolds number
    Re = rho * v * lc / mu

    return T,Tm,p,rho,ainf,M,mu,Pr,Kn,d,Re,R