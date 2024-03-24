# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 18:06:25 2023

@author: yanyan
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from math import acos, exp, sqrt, atan
from mpl_toolkits.mplot3d import Axes3D

def BEMT(TipLoss, N, Nb, sigma_function, theta_function, Cl_alpha, Cd_alpha):
    lambda_induced = np.zeros(N)
    delta_CT_t_ = np.zeros(N)
    delta_Cpi_t_ = np.zeros(N)
    delta_Cp0_t = np.zeros(N)
    delta_Cp_t = np.zeros(N)
    Cl = np.zeros(N)
    dr = 1/N

    for i in range(1, N+1):
        r_i = i * dr - dr / 2  # r_i
        F_i = 1
        term_inside_sqrt = 0  # 提前定义变量，确保在任何情况下都有值



        if sigma_function(r_i) * Cl_alpha != 0:
            term_inside_sqrt = 32 * F_i / (sigma_function(r_i) * Cl_alpha) * theta_function(r_i) * r_i
            if term_inside_sqrt >= 0:
                lambda_i = sigma_function(r_i) * Cl_alpha / (16 * F_i) * (sqrt(1 + term_inside_sqrt) - 1)
            else:
                lambda_i = 0  # 或者一些合适的默认值
        else:
            lambda_i = 0  # 或者一些合适的默认值        
        

        err = lambda_i
        n = 0

        while err >= 0.01 and n < 100:  # Added iteration limit to prevent infinite loop
            if TipLoss:
                F_i = F(Nb, r_i, lambda_i)
            else:
                F_i = 1
            
            if sigma_function(r_i) * Cl_alpha != 0:
                term = 32 * F_i / (sigma_function(r_i) * Cl_alpha) * theta_function(r_i) * r_i
                if term >= -1:
                    lambda_i_new = sigma_function(r_i) * Cl_alpha / (16 * F_i) * (sqrt(1 + term) - 1)
                    err = abs(lambda_i - lambda_i_new)
                    lambda_i = lambda_i_new
            n += 1

        alpha_i = theta_function(r_i) - atan(lambda_i / r_i) if r_i != 0 else 0
        Cl_i = Cl_alpha * alpha_i
        dCT = 1/2 * sigma_function(r_i) * Cl_i * r_i**2 * dr
        dCpi = lambda_i * dCT
        dCp0 = 1/2 * sigma_function(r_i) * Cd_alpha(alpha_i) * r_i**3 * dr

        lambda_induced[i-1] = lambda_i
        delta_CT_t_[i-1] = dCT / dr
        delta_Cpi_t_[i-1] = dCpi / dr
        delta_Cp0_t[i-1] = dCp0 / dr
        delta_Cp_t[i-1] = (dCpi + dCp0) / dr
        Cl[i-1] = Cl_i

    CT = sum(delta_CT_t_ * dr)
    CP = sum(delta_Cp_t * dr)

    return lambda_induced, delta_CT_t_, delta_Cpi_t_, delta_Cp0_t, delta_Cp_t, Cl, CT, CP


# Definition of the theta_0_solver function
def theta_0_solver(TipLoss, CT_req, Nb, theta_tw, sigma, Cl_alpha, Cd):
    def func_to_solve(theta0_j):
        theta = lambda r: theta0_j + theta_tw * r
        _, _, _, _, _, _, CT_j, _ = BEMT(TipLoss, 100, Nb, sigma, theta, Cl_alpha, Cd)
        return 6 * (CT_req - CT_j) / (sigma(0) * Cl_alpha)

    initial_guess = 6 * CT_req / (sigma(0) * Cl_alpha) - 3/4 * theta_tw + 3/2 * np.sqrt(CT_req / 2)
    theta0, = fsolve(func_to_solve, initial_guess)
    return theta0

# Definition of the F function
def F(Nb, r, lamda, epsilon=1e-6):
    factor = (Nb / 2) * (1 - r) / (lamda + epsilon)  # Add a small epsilon to avoid division by zero
    cos_factor = np.exp(-factor)
    # Ensure the value passed to arccos is within the valid range [-1, 1]
    cos_factor = np.clip(cos_factor, -1, 1)
    return (2 / np.pi) * np.arccos(cos_factor)
################################################################################
# Problem 5
# Constants and initial settings
Cl_alpha = 2 * np.pi
Cd_alpha = lambda alpha: 0.008
Nb = 2
theta_tw = -15  # degrees
taper_ratio = np.arange(1, 4)  # 1:1 to 3:1 taper ratios
sigma_weighted = 0.1
CT_req = 0.008
N = 80
r_BEMT = np.linspace(0.01, 1, N)
# Question 5 plots
plt.figure(10)
plt.suptitle(f'Question 5: #blades = {Nb}; C_T = {CT_req}') 

for gama in taper_ratio:
    sigma_initial = sigma_weighted / (1 + 3/4 * (1/gama - 1))
    sigma_iter = sigma_initial * (1/gama - 1)
    sigma_function = lambda r: sigma_initial + sigma_iter * r

    # Calculate theta_0 and set up theta_function
    theta_0 = theta_0_solver(True, CT_req, Nb, np.deg2rad(theta_tw), sigma_function, Cl_alpha, Cd_alpha)
    theta_function = lambda r: theta_0 + np.deg2rad(theta_tw) * r

    # BEMT calculations (use actual BEMT function)
    lambda_induced, delta_CT_t_, delta_Cpi_t_, delta_Cp0_t, delta_Cp_t, Cl, CT, CP = BEMT(True, N, Nb, sigma_function, theta_function, Cl_alpha, Cd_alpha)
    
    # Plotting for each taper ratio
    plt.subplot(2, 2, 1)
    plt.plot(r_BEMT, lambda_induced,'--', label=f'Taper Ratio = {gama}:1')
    plt.title("Inflow Distribution")
    plt.xlabel("r")
    plt.ylabel("λ")

    plt.subplot(2, 2, 2)
    plt.plot(r_BEMT, delta_CT_t_, '--', label=f'Taper Ratio = {gama}:1')
    plt.title("Thrust Distribution")
    plt.xlabel("r")
    plt.ylabel("dC_T/dr")

    plt.subplot(2, 2, 3)
    plt.plot(r_BEMT, delta_Cpi_t_,'--', label=f'Taper Ratio = {gama}:1')
    plt.title("Induced Torque Distribution")
    plt.xlabel("r")
    plt.ylabel("dC_{Q_i}/dr")

    plt.subplot(2, 2, 4)
    plt.plot(r_BEMT, Cl,'--', label=f'Taper Ratio = {gama}:1')
    plt.title("Lift Coeff. Distribution")
    plt.xlabel("r")
    plt.ylabel("C_L")

plt.legend()
plt.tight_layout()


# FM distribution calculations
taper_r = np.arange(0.5, 35.5, 0.5)  # Taper ratios
ang_tw = np.arange(-25, 6, 1)  # Twist angles in degrees
FM_fig = np.zeros((len(ang_tw), len(taper_r)))

for ang_index, twist_ang in enumerate(ang_tw):
    for b_index, gama_ in enumerate(taper_r):
        sigma_initial = sigma_weighted / (1 + 3/4 * (1/gama_ - 1))
        sigma_iter = sigma_initial * (1/gama_ - 1)
        sigma_function = lambda r: sigma_initial + sigma_iter * r

        theta_0 = theta_0_solver(True, CT_req, Nb, np.deg2rad(twist_ang), sigma_function, Cl_alpha, Cd_alpha)
        theta_function = lambda r: theta_0 + np.deg2rad(twist_ang) * r

        # BEMT calculations
        lambda_induced, delta_CT_t_, delta_Cpi_t_, Cp0_, delta_Cp_t, Cl, CT, CP = BEMT(True, N, Nb, sigma_function, theta_function, Cl_alpha, Cd_alpha)
        
        # Calculate FM and prevent imaginary results
        FM = np.real((CT ** (3/2) / sqrt(2)) / CP)
        FM_fig[ang_index, b_index] = FM

# Find maximum FM and corresponding twist and taper ratio
FM_max = np.max(FM_fig)
max_index = np.unravel_index(FM_fig.argmax(), FM_fig.shape)
twist_opt = ang_tw[max_index[0]]
taper_opt = taper_r[max_index[1]]

# Print optimal values
print(f"Maximum FM = {FM_max}")
print(f"Optimal angle of twist = {twist_opt} degrees")
print(f"Optimal Taper ratio = {taper_opt}")

# Plotting FM distribution
# Create meshgrids for ang_tw and taper_r
# ... [previous code]

TW, TR = np.meshgrid(ang_tw, taper_r, indexing='ij')

# Ensure that the shapes of TW, TR, and FM_fig are compatible
print("Shapes - TW:", TW.shape, "TR:", TR.shape, "FM_fig:", FM_fig.shape)

# Plotting FM distribution with corrected meshgrids
plt.figure()
ax = plt.axes(projection='3d')
ax.plot_surface(TW, TR, FM_fig, cmap='viridis', edgecolor='none')
ax.set_title('FM Distribution')
ax.set_xlabel('theta_{twist}')
ax.set_ylabel('Taper Ratio')
ax.set_zlabel('FM')
plt.show()
































