# -*- coding: utf-8 -*-
"""
Created on Fri Oct 10 15:17:09 2025

@author: alexa
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from scipy.integrate import quad
from scipy.optimize import curve_fit
from uncertainties import ufloat


font_size = 25


c = 299792.458


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
mu_err = df["mu_err"].values



# Part (ii) -------------------------------------------------------------------

# Convert from the distance modulus to the distance in Mega parsecs
dist_x = 10**(mu / 5.0 - 5.0)       

# Calculate errors
dist_err = np.log(10) / 5.0 * dist_x * mu_err



# Part (iii) ------------------------------------------------------------------

# Choose a mask value - 0.1 is suggested (only for part (iii))
# mask_var = 0.1 
# mask = df["redshift"] < mask_var


# z_mask  = df["redshift"][mask].values
# dist_mask = dist_x[mask]
# d_err_mask = dist_err[mask]


# coeffs, cov = np.polyfit(z_mask, dist_mask, 1, w=1/d_err_mask, cov=True)
# slope, intercept = coeffs
# slope_err, intercept_err = np.sqrt(np.diag(cov))

# # Convert slope to H0
# H0 = c / slope
# H0_err = c * slope_err / slope**2

# H_0_full = ufloat(H0, H0_err)

# print(f"H0 = {H_0_full} km/s/Mpc")


# plt.errorbar(z_mask, dist_mask, yerr=d_err_mask, fmt='.', alpha=0.6, label='Data (z<0.1)')
# plt.plot(z_mask, slope*z_mask, 'r', label=f'Fit: H0 = {H0:.1f} ± {H0_err:.1f} km/s/Mpc')
# plt.plot(z_mask, dist_mask, 'o', color='blue')
# plt.xlabel("Redshift (z)", fontsize = font_size)
# plt.ylabel("Distance (Mpc)", fontsize = font_size)
# plt.title("$Distance vs. redshift$", fontsize = font_size)
# plt.grid(True)

# plt.plot()




# Part (iv) -------------------------------------------------------------------
# '''
# We don't need to report \Omega_{Lambda,0} because of the identity that for flat space, 
# \Omega_{M} + \Omega_{Lambda,0} = 1. Hence, solving for one immediately gives the other.
# '''
# plt.figure()
# def d_fit(z, H0, Omega_m):
#     Omega_lambda = 1 - Omega_m
#     def integral(z_x):
#         E = np.sqrt(Omega_m * (1 + z_x)**3 + Omega_lambda)
#         return c / (H0 * E)
    
#     #z = np.atleast_1d(z)
#     D_C = np.array([quad(integral, 0, z_i)[0] for z_i in z])
#     return (1 + z) * D_C
    
# initial_vec = [65, 0.3]             
# bounds = ([10.0, 1.0], [200.0, 1.0])
# popt, pcov = curve_fit(d_fit, z, dist_x, p0=initial_vec, sigma=dist_err)

# H0_fit, Omega_m_fit = popt
# err = np.sqrt(np.diag(pcov))     
# H0_ferr, Omega_m_err = err

# H_full = ufloat(H0_fit, H0_ferr)
# omega_M_full = ufloat(Omega_m_fit, Omega_m_err)
# omega_Lambda_full = 1 - omega_M_full

# slope_2 = c / H0_fit
# print(f"H0 = {H_full} km/s/Mpc")
# print(f"Omega_m = {omega_M_full}")
# print(f"Omega_Lambda = {omega_Lambda_full} ")  


# plt.errorbar(z, dist_x, yerr=dist_err, fmt='.', alpha=0.6, label='Data (z<0.1)')
# plt.plot(z, slope_2*z, 'orange', label=f'Fit: H0_fit = {H0_fit:.1f} ± {H0_ferr:.1f} km/s/Mpc') 
# plt.xlabel("Redshift (z)", fontsize = font_size)
# plt.ylabel("Distance (Mpc)", fontsize = font_size)
# plt.title("Fitted Hubble constant for Distance vs redshift", fontsize = font_size)
# plt.legend()
# plt.plot()

# Part (v) --------------------------------------------------------------------
# If there was no dark energy, then Omega_lambda = 0


# plt.figure()
# def d_fit(z, H0, Omega_m):
#     def integral(z_x):
#         E = np.sqrt(Omega_m * (1 + z_x)**3)
#         return c / (H0 * E)
    
#     #z = np.atleast_1d(z)
#     D_C = np.array([quad(integral, 0, z_i)[0] for z_i in z])
#     return (1 + z) * D_C
    
# initial_vec = [70, 1]  
# bounds = ([10, 0.999], [200, 1])           

# popt, pcov = curve_fit(d_fit, z, dist_x, p0=initial_vec, sigma=dist_err, bounds=bounds)

# H0_fit, Omega_m_fit = popt
# err = np.sqrt(np.diag(pcov))     
# H0_ferr, Omega_m_err = err

# H_full = ufloat(H0_fit, H0_ferr)
# omega_M_full = ufloat(Omega_m_fit, Omega_m_err)
# omega_Lambda_full = 1 - omega_M_full

# slope_2 = c / H0_fit
# print(f"H0 = {H_full} km/s/Mpc")
# print(f"Omega_m = {omega_M_full}")
# print(f"Omega_Lambda = {omega_Lambda_full} ")  


# plt.errorbar(z, dist_x, yerr=dist_err, fmt='.', alpha=0.6, label='Data (z<0.1)')
# plt.plot(z, slope_2*z, 'orange', label=f'Fit: H0_fit = {H0_fit:.1f} ± {H0_ferr:.1f} km/s/Mpc') 
# plt.xlabel("Redshift (z)", fontsize = font_size)
# plt.ylabel("Distance (Mpc)", fontsize = font_size)
# plt.title("Fitted Hubble constant for Distance vs redshift", fontsize = font_size)
# plt.legend()
# plt.plot();

# Part (vi) -------------------------------------------------------------------


'''
From part (iv): H0 = 70.38+/-0.33 km/s/Mpc
                Omega_m = 0.297+/-0.019
                Omega_Lambda = 0.703+/-0.019 
'''
c = 299792.458
H0_x = 70.38
Omega_mx = 0.297
Omega_Lambdax = 0.703

def d_find(z):
    E = np.sqrt(Omega_mx * (1 + z)**3 + Omega_Lambdax)
    return c / (H0_x * E) 

z_max = max(z)
z_eval = np.linspace(0, z_max, 500)

d_phys = np.array([quad(d_find, 0, z_i)[0] for z_i in z_eval])

D_angular = d_phys / (1 + z_eval)

# plt.errorbar(z, dist_x, yerr=dist_err, fmt='.', alpha=0.6, label='Data')
# plt.plot(z_eval, D_angular, 'orange') 
# plt.show()


D_angular = D_angular[1:]

l_gal = 30 / 1000 
theta_rad = l_gal/ D_angular
theta_arcsec = theta_rad * (180/np.pi) * 3600

z_eval = z_eval[1:]


min_angle = np.min(theta_arcsec)
print(f"Minimum angular size: {min_angle:.3f} arcsec")

arsec_index = np.where(theta_arcsec == min_angle)
redshift_val = z_eval[arsec_index][0]
print(f"The redshift associated with this minimum angular size is {redshift_val}")

plt.plot(redshift_val, min_angle, 'b^', markersize=10, label='Most distant supernova')
plt.plot(z_eval, theta_arcsec, color='red', lw=3)
plt.xlabel('Redshift z',fontsize = font_size)
plt.ylabel('Angular size [arcsec]',fontsize = font_size)
plt.title('Minimum Angular Size of a 30 kpc Galaxy',fontsize = font_size)
plt.grid(True)
plt.legend()
plt.show()





def lookback(t):
    E = np.sqrt(Omega_mx * (1 + t)**3 + Omega_Lambdax)
    return 1 / (H0_x * (1 + t) * E) 
    
t_lookback, _ = quad(lookback, 0, redshift_val)

conversion = 3.08*10**19 * 3.17*10**-17
print(f"The lookback time is therefore {t_lookback * conversion} Gy")

print(""" This is a very large lookback time, and occured when the Universe was just 4.45 billions years old. 
At these high redshifts, galaxies could be more luminous than the supernovae making it difficult to see 
these supernovae in reality.""")






