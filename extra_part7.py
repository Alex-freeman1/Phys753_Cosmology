# -*- coding: utf-8 -*-
"""
Created on Thu Oct 16 20:04:47 2025

@author: alexa
"""


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from scipy.integrate import quad
from scipy.optimize import curve_fit
from uncertainties import ufloat

# Part (vii - optional)
c = 299792.458
font_size = 25

# # Load data

C_mu = np.loadtxt("https://supernova.lbl.gov/Union/figures/SCPUnion2.1_covmat_sys.txt")
url = "https://supernova.lbl.gov/Union/figures/SCPUnion2.1_mu_vs_z.txt"

# Load the table
df = pd.read_csv(url, comment='#', sep=r'\s+', header=None)
df.columns = [
    "supernova_id",
    "redshift",
    "mu",
    "mu_err",
    "low_mass_prob"
]


z = df["redshift"].values
mu = df["mu"].values

     
# Convert modulus to d_L
dist_x = 10**(mu / 5.0 - 5.0)  

# Jacobian diagonal
J = (np.log(10) / 5) * dist_x

# Propagate covariance
C_dL = np.outer(J, J) * C_mu
dist_err = np.sqrt(np.diag(C_dL))



# Choose a mask value - 0.1 is suggested (only for part (iii))
mask_var = 0.1
mask = z < mask_var

z_mask = z[mask]
dist_mask = dist_x[mask]
d_err_mask = dist_err[mask]


coeffs, cov = np.polyfit(z_mask, dist_mask, 1, w=1/d_err_mask, cov=True)
slope, intercept = coeffs
slope_err, intercept_err = np.sqrt(np.diag(cov))

# Convert slope to H0
H0 = c / slope
H0_err = c * slope_err / slope**2

H_0_full = ufloat(H0, H0_err)

print(f"H0 = {H_0_full} km/s/Mpc")


plt.errorbar(z_mask, dist_mask, yerr=d_err_mask, fmt='.', alpha=0.6, label='Data (z<0.1)')
plt.plot(z_mask, slope*z_mask, 'r', label=f'Fit: H0 = {H0:.1f} ± {H0_err:.1f} km/s/Mpc')
plt.plot(z_mask, dist_mask, 'o', color='blue')
plt.xlabel("Redshift (z)", fontsize = font_size)
plt.ylabel("Distance (Mpc)", fontsize = font_size)
plt.title("$Distance vs. redshift$", fontsize = font_size)
plt.grid(True)

plt.plot()

plt.figure()
def d_fit(z, H0, Omega_m):
    Omega_lambda = 1 - Omega_m
    def integral(z_x):
        E = np.sqrt(Omega_m * (1 + z_x)**3 + Omega_lambda)
        return c / (H0 * E)
    
    #z = np.atleast_1d(z)
    D_C = np.array([quad(integral, 0, z_i)[0] for z_i in z])
    return (1 + z) * D_C
    
initial_vec = [65, 0.3]             
bounds = ([10.0, 1.0], [200.0, 1.0])
popt, pcov = curve_fit(d_fit, z, dist_x, p0=initial_vec, sigma=dist_err)

H0_fit, Omega_m_fit = popt
err = np.sqrt(np.diag(pcov))     
H0_ferr, Omega_m_err = err

H_full = ufloat(H0_fit, H0_ferr)
omega_M_full = ufloat(Omega_m_fit, Omega_m_err)
omega_Lambda_full = 1 - omega_M_full

slope_2 = c / H0_fit
print(f"H0 = {H_full} km/s/Mpc")
print(f"Omega_m = {omega_M_full}")
print(f"Omega_Lambda = {omega_Lambda_full} ")  


plt.errorbar(z, dist_x, yerr=dist_err, fmt='.', alpha=0.6, label='Data (z<0.1)')
plt.plot(z, slope_2*z, 'orange', label=f'Fit: H0_fit = {H0_fit:.1f} ± {H0_ferr:.1f} km/s/Mpc') 
plt.xlabel("Redshift (z)", fontsize = font_size)
plt.ylabel("Distance (Mpc)", fontsize = font_size)
plt.title("Fitted Hubble constant for Distance vs redshift", fontsize = font_size)
plt.legend()
plt.plot();


# Part (viiI - optional)



