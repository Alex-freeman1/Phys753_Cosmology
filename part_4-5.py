# -*- coding: utf-8 -*-
"""
Created on Thu Oct 16 13:12:25 2025

@author: alexa
"""

import numpy as np

from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


font_size = 25


# Constants {speed of light (km/s), hubble's constant, matter density, and dark energy density}
c = 299792.458  
H_0 = 70  
Omega_m = 0.3
Omega_lambda = 0.7



# Part 4 ------------------------------------------------------------------------------------------------------------

# Define luminosity distance
def distance_c(z, D):
    E = np.sqrt(Omega_m * (1 + z)**3 + Omega_lambda)
    return c / (H_0 * E)

# Declare the limits of the integral - max redshift observed
z_max = 8
z_eval = np.linspace(0, z_max, 500)
z_span = (z_eval[0], z_eval[-1])
sol = solve_ivp(distance_c, z_span, [0], t_eval=z_eval, method='RK45', rtol=1e-10, atol=1e-10)

# Compute luminosity distance
D_L = (1 + z_eval) * sol.y[0]


# Plotting the luminosity distance using the \Lambda CDM model
plt.plot(z_eval, D_L, label=" $(\Omega_{M}, \Omega_{\Lambda}) = (0.3,0.7). \Lambda CDM$", color ='blue', linewidth=4)

plt.title(r'Luminsoity distance vs redshift', fontsize = font_size)
plt.xlabel("Redshift $z$", fontsize = font_size)
plt.ylabel(r"Distance $D_L$ [Mpc]", fontsize = font_size)
plt.grid(True)
plt.legend(fontsize = font_size)
plt.show()



# Part 5 ------------------------------------------------------------------------------------------------------------
plt.figure() 

# Using exact equations found in part (3)
DL_denergy = c / H_0 * z_eval * (1 + z_eval)
DL_matter = 2 * c / H_0 * (1 + z_eval) * (1 - 1 / (1 + z_eval)**0.5)


# Plotting all three curves (Same D_L is used as in part (4))
plt.plot(z_eval, DL_denergy, label=" $(\Omega_{M}, \Omega_{\Lambda}) = (0,1). Dark energy-dominated$", color ='orange', linewidth=4)
plt.plot(z_eval, D_L, label=" $(\Omega_{M}, \Omega_{\Lambda}) = (0.3,0.7). \Lambda CDM$", color ='blue', linewidth=4)
plt.plot(z_eval, DL_matter, label=" $(\Omega_{M}, \Omega_{\Lambda}) = (1,0). Matter-dominated$", color ='red', linewidth=4)

plt.title(r'Luminsoity distance calculations vs redshift', fontsize = font_size)
plt.xlabel("Redshift $z$", fontsize = font_size)
plt.ylabel(r"Distance $D_L$ [Mpc]", fontsize = font_size)
plt.grid(True)
plt.legend(fontsize = font_size)
plt.show()





    
