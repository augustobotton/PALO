import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

def res3body(t,X):
    m = 0.01215  # massa não dimensional do segundo primário
    x,y,xp1,xp2 = X  

    # distâncias dos primários
    r1 = np.sqrt((x - m) ** 2 + y ** 2)
    r2 = np.sqrt((x + 1 - m) ** 2 + y ** 2)

    # Equações de movimento:
    xp3 = x + 2 * xp2 - (1 - m) * (x - m) / r1 ** 3 - m * (x + 1 - m) / r2 ** 3
    xp4 = y - 2 * xp1 - (1 - m) * y / r1 ** 3 - m * y / r2 ** 3

    return [xp1,xp2,xp3,xp4]

# vetor de tempo
t_span = (0, 8)

# condições iniciais
comp_velocidade = [
    [0.1, 0, 0, 0.5],     # (a)
    [0.1, 0, -4, 1],      # (b)
    [0.1, 0, -3.35, 3],   # (c)
    [0.1, 0, -3.37, 3],   # (d)
    [0.1, 0, -3.4, 3],    # (e)
    [0.1, 0, -3.5, 3],    # (f)
    [0.1, 0, -3.6, 3],    # (g)
]


#X0 = np.array([0.1, 0, -3.5, 3])
for i, X0 in enumerate(comp_velocidade, 1):
    sol = solve_ivp(res3body, t_span, X0, t_eval=np.linspace(t_span[0], t_span[1], 1000000),dense_output=True,method='RK45')
    x, y = sol.y[0], sol.y[1]
    
    plt.figure()  
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title(f"Trajetória da Espaçonave: Caso {chr(96+i)}")
    plt.plot(x, y, label='Trajetória')
    #plt.plot(-0.99,0,'x',label='Lua')
    #plt.plot(0.0125,0,'x',label='Terra')
    plt.legend()
    
plt.show()