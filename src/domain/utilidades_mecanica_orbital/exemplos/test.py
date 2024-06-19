import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation

# Gravitational constant
G = 6.67430e-11  # in m^3 kg^-1 s^-1

# Define the system of ODEs for the n-body problem
def n_body_equations(w, t, masses):
    n = len(masses)
    positions = w[:n * 3].reshape((n, 3))
    velocities = w[n * 3:].reshape((n, 3))
    derivatives = np.zeros_like(w)

    for i in range(n):
        # Calculate accelerations due to gravity
        acceleration = np.zeros(3)
        for j in range(n):
            if i != j:
                r = positions[j] - positions[i]
                acceleration += G * masses[j] * r / np.linalg.norm(r) ** 3

        # Update velocity derivatives
        derivatives[i * 3:(i + 1) * 3] = velocities[i]

        # Update position derivatives
        derivatives[n * 3 + i * 3:n * 3 + (i + 1) * 3] = acceleration

    return derivatives

# Initialize parameters
n_bodies = 5  # Number of bodies (Sun, Mercury, Venus, Earth, Mars)
masses = np.array([1.989e30, 3.285e23, 4.867e24, 5.972e24, 6.39e23])  # Masses in kg

# Initial positions in meters (approximate)
positions = np.array([
    [0, 0, 0],  # Sun
    [5.79e10, 0, 0],  # Mercury
    [1.082e11, 0, 0],  # Venus
    [1.496e11, 0, 0],  # Earth
    [2.279e11, 0, 0]  # Mars
])

# Initial velocities in meters/second (approximate)
velocities = np.array([
    [0, 0, 0],  # Sun
    [0, 47360, 0],  # Mercury
    [0, 35020, 0],  # Venus
    [0, 29780, 0],  # Earth
    [0, 24130, 0]  # Mars
])

# Flatten the initial conditions into a single array
initial_conditions = np.hstack((positions.flatten(), velocities.flatten()))

# Time span
t = np.linspace(0, 6.154e7, 1000)  # One year in seconds, with 1000 time points

# Solve the ODEs
solution = odeint(n_body_equations, initial_conditions, t, args=(masses,))

# Extract positions from the solution
positions_sol = solution[:, :n_bodies * 3].reshape((len(t), n_bodies, 3))

# Create the figure and axis for the animation
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Set plot limits (optional)
ax.set_xlim([-2e11, 2e11])
ax.set_ylim([-2e11, 2e11])
ax.set_zlim([-2e11, 2e11])

# Initialize the markers for each body
markers = [ax.plot([], [], [], marker='o')[0] for _ in range(n_bodies)]

# Variables to control the animation
paused = False
interval = 50  # Interval between frames in milliseconds

# Update function for the animation
def update(num):
    for i in range(n_bodies):
        markers[i].set_data(positions_sol[num, i, 0], positions_sol[num, i, 1])
        markers[i].set_3d_properties(positions_sol[num, i, 2])
    return markers

# Pause/play function
def toggle_pause(event):
    global paused
    if event.key == ' ':
        paused = not paused
        if paused:
            ani.event_source.stop()
        else:
            ani.event_source.start()

# Speed control function
def adjust_speed(event):
    global interval
    if event.key == 'up':
        interval = max(10, interval - 10)  # Increase speed (decrease interval)
        ani.event_source.interval = interval
    elif event.key == 'down':
        interval += 10  # Decrease speed (increase interval)
        ani.event_source.interval = interval

# Create the animation
ani = FuncAnimation(fig, update, frames=len(t), interval=interval, blit=True)

# Connect the event handlers
fig.canvas.mpl_connect('key_press_event', toggle_pause)
fig.canvas.mpl_connect('key_press_event', adjust_speed)

# Display the animation
plt.show()

