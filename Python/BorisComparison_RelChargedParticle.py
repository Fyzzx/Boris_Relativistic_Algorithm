# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 14:47:38 2023

@author: pry81
"""

import numpy as np
import matplotlib.pyplot as plt
for Boris in range(1,4):
    # Constants
    q = 1.0  # charge
    m = 1.0  # mass
    c = 1.0  # speed of light
    
    # Initial conditions
    init_pos = np.array([0.9, 0.0, 0.0])  # Initial position (at the +1/2 position)
    init_u = np.array([0.1, 0.0, 0.0])    # Initial velocity multiplied by gamma
    
    # Time parameters
    dt = np.pi / 1E2      # Time step
    total_time = 1E3      # Total time to run simulation
    num_steps = int(np.ceil(total_time / dt))  # Number of time steps
    
    # Preallocate arrays
    pos = np.zeros((num_steps+1, 3))
    u = np.zeros((num_steps+1, 3))
    e_Field = np.zeros((num_steps+1, 3))
    b_Field = np.zeros((num_steps+1, 3))
    
    # Initialize variables
    pos[0, :] = init_pos  # Position at the 1/2 steps
    u[0, :] = init_u      # Velocity times gamma at nn = 1
    
    # Iterative loop over "nn" time steps
    for nn in range(num_steps):
    
        # Find the current values for E and B fields at the n+1/2 value
        pos_plusHalf = pos[nn, :]
        
        # Find E field (nn+1/2)
        eFactor = (0.01 * (pos_plusHalf[0]**2 + pos_plusHalf[1]**2)**(-3/2))
        e_Field[nn, :] = eFactor * np.array([pos_plusHalf[0], pos_plusHalf[1], 0])
        
        # Find B field (nn+1/2)
        b_Field[nn, :] = np.array([0, 0, np.sqrt(pos_plusHalf[0]**2 + pos_plusHalf[1]**2)])
        
        # 1st half step push of particle due to E field [Equation 3]
        u_minus = u[nn, :] + (q/m) * (dt/2) * e_Field[nn, :]
    
        # Rotation piece - Boris (2 step rotation) [Equation 4] u_minus -> u_plus
        u_minus_mag = np.linalg.norm(u_minus)
        gamma_minus = np.sqrt(1 + (u_minus_mag/c)**2)
        B_mag = np.sqrt(b_Field[nn, 0]**2 + b_Field[nn, 1]**2 + b_Field[nn, 2]**2)
        theta = (q * dt * B_mag) / (m * gamma_minus)
    
        if Boris < 3:
            if Boris == 1:  # Boris A
                t_vector = np.tan(theta / 2) * (b_Field[nn, :] / B_mag)
            else:  # Boris B
                t_vector = (theta / 2) * (b_Field[nn, :] / B_mag)
    
            # Step 1 of rotation u_minus -> u_prime [Equation 8]
            u_prime = u_minus + np.cross(u_minus, t_vector)
            t_squared = np.dot(t_vector, t_vector)  # Needed for Equation 9
    
            # Updating the final piece of the rotation u_prime -> u_plus [Equation 9]
            u_plus = u_minus + (2 / (1 + t_squared)) * np.cross(u_prime, t_vector)
        else:  # Boris C
            # Equation 11
            u_parallel_minus = np.dot(u_minus, (b_Field[nn, :] / B_mag)) * (b_Field[nn, :] / B_mag)
            
            # Equation 12
            u_plus = u_parallel_minus + (u_minus - u_parallel_minus) * np.cos(theta) + \
                     np.cross(u_minus, (b_Field[nn, :] / B_mag)) * np.sin(theta)
    
        # Final E field half time step push [Equation 5]
        u[nn + 1, :] = u_plus + (q * dt / (2 * m)) * e_Field[nn, :]
    
        # Update the position [Equation 1]
        u_mag = np.linalg.norm(u[nn + 1, :])
        gamma = np.sqrt(1 + (u_mag/c)**2)
        pos[nn + 1, :] = dt * u[nn + 1, :] / gamma + pos[nn, :]
        
        if Boris == 1:
            posA = pos
        
        if Boris == 2:
            posB = pos

        if Boris == 3:
            posC = pos
        
# Visualization (similar to MATLAB)
FS = 25

plt.figure(456)
posDiffAC = posA - posC
plt.plot(posDiffAC[:, 0], posDiffAC[:, 1])
plt.xlabel('X Position')
plt.ylabel('Y Position')
plt.title('Position difference between Boris A and C', fontsize=FS)
plt.plot(0, 0, 'k.', markersize=30)
ax = plt.gca()
ax.tick_params(axis='both', which='major', labelsize=FS)

plt.figure(789)
posDiffBC = posB - posC
plt.plot(posDiffBC[:, 0], posDiffBC[:, 1])
plt.xlabel('X Position')
plt.ylabel('Y Position')
plt.title('Position difference between Boris B and C', fontsize=FS)
ax = plt.gca()
ax.tick_params(axis='both', which='major', labelsize=FS)

plt.figure(123)
plt.plot(posA[:, 0], posA[:, 1], 'k')
plt.plot(posB[:, 0], posB[:, 1], 'r')
plt.plot(posC[:, 0], posC[:, 1], 'g')
plt.xlabel('X Position')
plt.ylabel('Y Position')
plt.axis('equal')
plt.title('Relativistic motion of charged particle using Boris A, B, and C Methods', fontsize=FS)
plt.legend(['Boris A', 'Boris B', 'Boris C'], fontsize=FS)
ax = plt.gca()
ax.tick_params(axis='both', which='major', labelsize=FS)

plt.show()
