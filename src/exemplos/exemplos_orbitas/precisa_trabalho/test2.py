import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation
from jplephem.spk import SPK
from datetime import datetime, timedelta

# Load the JPL ephemeris data
kernel = SPK.open('C:/Users/gt_po/Documents/tcc/mecvooespacial/src/de440t.bsp')

# Define the simulation time span
start_date = datetime(2024, 1, 1)
end_date = datetime(2025, 1, 1)
num_steps = 1000
times = np.linspace(0, (end_date - start_date).total_seconds(), num_steps)


# Define a function to get position and velocity data from the ephemeris
def get_ephemeris_data(body_id, time_seconds):
	time_jd = 2451545.0 + time_seconds / 86400.0  # Convert to Julian date
	position, velocity = kernel[0, body_id].compute_and_differentiate(time_jd)
	return position, velocity


# Define the bodies and their IDs in the ephemeris file
bodies = {
	'Sun': 10,
	'Mercury': 1,
	'Venus': 2,
	'Earth': 3,
	'Mars': 4
}
n_bodies = len(bodies)

# Get initial conditions from the ephemeris
initial_positions = []
initial_velocities = []
for body in bodies.values():
	position, velocity = get_ephemeris_data(body, 0)
	initial_positions.append(position)
	initial_velocities.append(velocity)

# Convert to numpy arrays
initial_positions = np.array(initial_positions)
initial_velocities = np.array(initial_velocities)

# Create the figure and axis for the animation
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Set plot limits (optional)
ax.set_xlim([-2e11, 2e11])
ax.set_ylim([-2e11, 2e11])
ax.set_zlim([-2e11, 2e11])

# Initialize the lines for each body and set labels for the legend
lines = [ax.plot([], [], [], marker='o', label=body)[0] for body in bodies.keys()]

# Add the legend to the plot
ax.legend()


# Update function for the animation
def update(num):
	time_seconds = times[num]
	positions = []
	for body in bodies.values():
		position, _ = get_ephemeris_data(body, time_seconds)
		positions.append(position)

	positions = np.array(positions)

	for i, line in enumerate(lines):
		line.set_data(positions[i, 0], positions[i, 1])
		line.set_3d_properties(positions[i, 2])
	return lines


# Create the animation
ani = FuncAnimation(fig, update, frames=num_steps, interval=50, blit=True)

# Display the animation
plt.show()
