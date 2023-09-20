# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import matplotlib.pyplot as plt

# Constants
q = 1.0  # Particle charge
m = 1.0  # Electron mass
c = 1.0  # Speed of light

# Initial conditions
init_pos = np.array([0.9, 0.0, 0.0])  # Initial position (at the +1/2 position)
init_u = np.array([0.1, 0.0, 0.0])  # Initial velocity multiplied by gamma

# Time parameters
dt = np.pi / 1E2  # Time step
total_time = 1.0E3  # Total time to run simulation
num_steps = int(np.ceil(total_time / dt))  # Number of time steps

# Preallocate arrays
pos = np.zeros((num_steps+1, 3))
u = np.zeros((num_steps+1, 3))
e_Field = np.zeros((num_steps+1, 3))
b_Field = np.zeros((num_steps+1, 3))

# Initialize variables
pos[0, :] = init_pos  # Position at the 1/2 steps
u[0, :] = init_u  # Velocity times gamma at nn = 1

# Iterative loop over "nn" time steps
for nn in range(num_steps):

    # Find the current values for E and B fields at the n+1/2 value
    pos_plusHalf = pos[nn, :]
    # Find E field (nn+1/2)
    eFactor = (0.01 * (pos_plusHalf[0]**2 + pos_plusHalf[1]**2)**(-3/2))
    e_Field[nn, :] = eFactor * np.array([pos_plusHalf[0], pos_plusHalf[1], 0])
    # Find B field (nn+1/2)
    b_Field[nn, :] = [0, 0, np.sqrt(pos_plusHalf[0]**2 + pos_plusHalf[1]**2)]

    # 1st half step push of particle due to E field [Equation 3]
    u_minus = u[nn, :] + (q/m) * (dt/2) * e_Field[nn, :]

    # Rotation piece - Boris (2 step rotation) [Equation 4] u_minus -> u_plus
    # Boris A = [Eqs. 6, 7a, 8, 9]

    # Compute magnitude of u_minus [needed for gamma_minus calculation]
    u_minus_mag = np.sqrt(np.sum(u_minus**2))

    # Compute gamma at the u_minus 'velocity' [inline equation between eqs. 2 & 3]
    gamma_minus = np.sqrt(1 + (u_minus_mag/c)**2)

    # Magnitude of B field at n+1/2
    B_mag = np.sqrt(np.sum(b_Field[nn, :]**2))

    # [Equation 6] phase angle
    theta = (q * dt * B_mag) / (m * gamma_minus)

    # [Equation 7a]
    t_vector = np.tan(theta / 2) * (b_Field[nn, :] / B_mag)

    # [Equation 8] step 1 of rotation u_minus -> u_prime
    u_prime = u_minus + np.cross(u_minus, t_vector)

    # Needed for Equation 9
    t_squared = np.dot(t_vector, t_vector)

    # [Equation 9] updating the final piece of the rotation u_prime -> u_plus
    u_plus = u_minus + (2 / (1 + t_squared)) * np.cross(u_prime, t_vector)

    # Final E field half time step push [Equation 5]
    u[nn + 1, :] = u_plus + (q * dt / (2 * m)) * e_Field[nn, :]

    # Update the position [Equation 1]

    # Needed to compute gamma for equation 1
    u_mag = np.linalg.norm(u[nn + 1, :])
    gamma = np.sqrt(1 + (u_mag/c)**2)

    # [Equation 1]
    pos[nn + 1, :] = dt * u[nn + 1, :] / gamma + pos[nn, :]

# Plot the results
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(pos[:, 0], pos[:, 1], pos[:, 2])
ax.set_xlabel('X Position')
ax.set_ylabel('Y Position')
ax.set_zlabel('Z Position')
ax.set_title(r'Relativistic motion of electron using Boris A Method, $\Delta t=$' + str(round(dt,6)) + 's')
plt.show()
