
'''

1. A rotor has a rectangular blade planform, a solidity of 0.1, ideal blade twist, and operates in hover at a thrust coefficient of 0.008. 
By numerical means, using the combined blade element momentum theory (BEMT), compare the predicted inflow ratio across the span of the blade with the exact (analytic) result (uniform inflow result). 
Plot continuous line for exact result and points (markers) for numerical result. Show also the radial distributions of thrust, induced torque, lift coefficient, and the compare the solutions obtained using the exact theory. 
Compare the numerical calculations of the integrated thrust and induced power with exact results for this rotor over the range 0 < CT < 0.01.

2. Using BEMT, implement Prandtl’s tip-loss function. For a rotor of solidity 0.1 and with no blade twist, show the effects of the tip-loss effect on the induced inflow, thrust distribution,
lift coefficient distribution, and torque distribution across the blade span. Neglect profile drag in this case. 
Compare the results for rotors with 2 and 4 blades at the same thrust coefficient of 0.008. 
Also, calculate the induced power factor (κ) for the two rotors, with and without the tip-loss effect and comment on your results.

3. Using BEMT, show the effect of increasing linear twist on the variations in inflow, thrust, induced power, profile power and lift 
coefficient across the span of a rotor with 4 blades of rectangular planform and solidity 0.1, and operating at a thrust coefficient of 0.008. 
Assume Cd = 0.011. Include Prandtl’s tip-loss effects. All of the results should be compared at a constant thrust (disk loading).

4. For the same rotor parameters as used in Question 3, compute the power required to hover versus thrust coefficient for different values of linear blade twist. 
Also compute the figure of merit (FM) and induced power factor (κ) versus thrust coefficient. 
Assume that the profile drag is given by Cd = 0.011 – 0.025α + 0.65α2. Show your results graphically.
Include the Prandtl’s tip-loss function in the calculations. Comment on your results.

5. Tapering the blade planform can have a powerful impact on the profile power of the rotor.
Using BEMT, show the effect of increasing taper on the variations in inflow, thrust, power and lift coefficient across the span of a rotor with 2 blades of linearly tapered planform and with a thrust weighted solidity of 0.1. 
The blades can be assumed to have 15° of linear nose-down (negative) twist. Consider blade taper ratios of 1:1, 2:1, and 3:1. 
Include Prandtl’s tip-loss effects. Remember that all of the results should be compared at a constant thrust (disk loading) and thrust weighted solidity. 
Determine an optimum combination of twist and taper that will give a maximum figure of merit (show your results using a contour plot or a surface plot).

'''
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from math import sqrt, atan, pi, radians
# from mpl_toolkits.mplot3d import Axes3D


def BEMT(TipLoss, N, Nb, sigma_function, theta_function, Cl_alpha, Cd_alpha):
    '''
    定义参数包括:
        TipLoss: 尖端损失效应
        N: 分割数量
        Nb: 叶片数量
        sigma_function: 实度函数
        theta_function: 扭转角函数
        Cl_alpha: 升力系数
        Cd_alpha: 阻力系数
    '''
    
    # initialization
    lambda_induced = np.zeros(N)
    delta_CT_t_ = np.zeros(N)
    delta_Cpi_t_ = np.zeros(N)
    delta_Cp0_t = np.zeros(N)
    delta_Cp_t = np.zeros(N)
    Cl = np.zeros(N)
    dr = 1/N    # the width of each segment


    for i in range(1, N+1): # Loop through each segment, 对每个分段进行循环
        
        r_i = i * dr - dr / 2  # r_i: 当前分段的半径位置
        F_i = 1                # correction factor: 初始化修正因子
        term_inside_sqrt = 0   # 提前定义变量，确保在任何情况下都有值,初始化平方根内部的项

        if sigma_function(r_i) * Cl_alpha != 0:
        # 如果实度函数和升力斜率乘积不为零:
            
            term_inside_sqrt = 32 * F_i / (sigma_function(r_i) * Cl_alpha) * theta_function(r_i) * r_i
            
            if term_inside_sqrt >= 0:
                
                lambda_i = sigma_function(r_i) * Cl_alpha / (16 * F_i) * (sqrt(1 + term_inside_sqrt) - 1)
                # inflow equation 诱导流量
                
            else:
                lambda_i = 0  # 或者一些合适的默认值
        
        else:
            lambda_i = 0  # 或者一些合适的默认值        
        

        err = lambda_i   # 初始化误差
        n = 0            # 初始化迭代次数

        while err >= 0.01 and n < 100:  # Added iteration limit to prevent infinite loop
            if TipLoss: # 尖端损失效应
                F_i = F(Nb, r_i, lambda_i) # Update F
            else:
                F_i = 1
            
            if sigma_function(r_i) * Cl_alpha != 0: 
                term = 32 * F_i / (sigma_function(r_i) * Cl_alpha) * theta_function(r_i) * r_i
                
                if term >= -1:
                    lambda_i_new = sigma_function(r_i) * Cl_alpha / (16 * F_i) * (sqrt(1 + term) - 1)
                    err = abs(lambda_i - lambda_i_new)
                    lambda_i = lambda_i_new
            n += 1    # 迭代次数加一

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
    
    initial_guess = 6 * CT_req / (sigma(0) * Cl_alpha) - 3/4 * theta_tw + 3/2 * sqrt(CT_req / 2)
    theta0, = fsolve(func_to_solve, initial_guess)
    return theta0


def F(Nb, r, lamda, epsilon=1e-6):
    factor = (Nb / 2) * (1 - r) / (lamda + epsilon)  # Add a small epsilon to avoid division by zero
    cos_factor = np.exp(-factor)
    # Ensure the value passed to arccos is within the valid range [-1, 1]
    cos_factor = np.clip(cos_factor, -1, 1)
    return (2 / np.pi) * np.arccos(cos_factor)


#######################################################################################################
# Problem 1:
# Constants
Nb = 2
Cl_alpha = 2 * np.pi
Cd = lambda alpha: 0  # assuming Cd is a function of alpha but returns 0 in this case
CT = 0.008  # Thrust coefficient, make sure this is defined before using it in lambda_analy
sigma = lambda r: 0.1  # solidity function

# Ideal Twist Calculation
theta_tip = 4 * CT / (0.1 * Cl_alpha) + np.sqrt(CT / 2)
theta = lambda r: theta_tip / r  # twist function

# Analytical Approach:
r_analy = np.arange(0.02, 1.02, 0.02)  # start from 0.02 to avoid division by zero
lambda_analy = sqrt(CT / 2) + 0 * r_analy
Cl_analy = 4 * CT / 0.1 / r_analy
Cl_analy[0] = 62  # to handle division by zero at r=0
CT_analy = 2 * CT * r_analy
Cp_induced_analy = lambda_analy * CT_analy
# Cp_induced_analy = np.sqrt(CT / 2) * CT_analy

# The CT_list_analy and Cp_induced_list_analy for the last subplot
CT_list_analy = np.arange(0, 0.0102, 0.0002)
# Cp_induced_list_analy = CT_list_analy ** (3 / 2) / np.sqrt(2) rdyjgukidtyik
Cp_induced_list_analy = np.power(CT_list_analy, 1.5) / sqrt(2)
#######################################################################################################

N = 50
r_BEMT = np.linspace(0.01, 1, N)
lambda_induced_list, delta_CT_t, delta_Cpi_t, _, _, Cl_list, _, CP_list = BEMT(
    False, N, Nb, sigma, theta, Cl_alpha, Cd)


#######################################################################################################

def plot_subplot(ax, x1, y1, x2, y2, title, xlabel, ylabel, legend_labels):
    ax.plot(x1, y1, 'C0', label=legend_labels[0], linewidth=1.2)  # 'k' is for black line
    ax.plot(x2, y2, 'o', label=legend_labels[1], color='orange', linestyle='')    
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.legend()
    
    
# Create figure and axes objects for subplots with a 2x3 grid
fig, axs = plt.subplots(2, 3, figsize=(18, 12))
fig.suptitle('Problem 1:')


# Plotting on the 2x3 grid
#plot 1
plot_subplot(axs[0, 0], r_analy, lambda_analy, r_BEMT, lambda_induced_list, 
    "f(r) = λ \n Radial distribution of inflow for the BEMT test case", 
    "Non-dimensional radial position, r", 
    "Induced inflow coefficient, λ", 
    ['Theory (Exact solution)', 'BEMT(Numerical solution)'])
axs[0, 0].set_ylim([0, 0.2])

# plot 2 前面的是exact，后面的是BEMT
plot_subplot(axs[0, 1],  r_analy, CT_analy,r_BEMT, delta_CT_t, 
    "$f(r) = C_{T}$ \n Radial distribution of thrust for the BEMT test case", 
    "Non-dimensional radial position, r",
    "$ dC_{T}/dr $", 
    ['Theory (Exact solution)', 'BEMT(Numerical solution)'])

# plot 3
plot_subplot(axs[0, 2], r_analy, Cp_induced_analy, r_BEMT, delta_Cpi_t, 
             "f(r) = ${C_Q}$ \n Radial distribution of induced power", 
             "Non-dimensional radial position, r",
             "$ dC_{Q}/dr $ ", 
             ['Theory (Exact solution)', 'BEMT(Numerical solution)'])

# Plot 4
plot_subplot(axs[1, 0], r_analy, Cl_analy, r_BEMT, Cl_list, 
             "f(r) = ${C_L}$ \n Radial distribution of lift coefficient for the BEMT test case", 
             "Non-dimensional radial position, r",
             "Lift coefficient, $ C_{l} $", 
             ['Theory (Exact solution)', 'BEMT(Numerical solution)'])

# Plot 5: A more fitted version of the lift coefficient distribution
plot_subplot(axs[1, 1], r_analy, Cl_analy, r_BEMT, Cl_list, 
             "f(r) = ${C_L}$ \n Refined Analysis (Enhanced Fit Version) \n Radial distribution of lift coefficient for the BEMT test case", 
             "Non-dimensional radial position, r",
             "Lift coefficient, $ C_{l} $", 
             ['Theory (Exact solution)', 'BEMT(Numerical solution)'])
axs[1, 1].set_ylim([0, 2])


# Predicted induced power versus thrust compared to analytic solution
# Pre-calculate values for the sixth plot 不能移动位置

# CT_list_analy = np.arange(0, 0.0102, 0.0002)
# Cp_induced_list_analy = CT_list_analy ** (3 / 2) / np.sqrt(2)
CT_val_ = np.arange(0.001, 0.0101, 0.001)
Cpi_val_ = []

for CT in CT_val_:
    theta_tip = 4 * CT / (sigma(0) * Cl_alpha) + np.sqrt(CT / 2)  # ideal twist => Solve for θ_tip
    theta = lambda r: theta_tip / r
    _, _, delta_Cpi_t, _, _, _, _, CP_list = BEMT(False, N, Nb, sigma, theta, Cl_alpha, Cd)
    Cpi_val_.append(sum(delta_Cpi_t / N))

# Convert list to numpy array for consistent data handling
Cpi_val_ = np.array(Cpi_val_)


# Plot 6
plot_subplot( axs[1, 2],   CT_list_analy, Cp_induced_list_analy,  CT_val_,  Cpi_val_, 
    "$f(C_{T}) = C_{P_i}$ \n Predicted induced power versus thrust compared to analytic solution", 
    "$C_{T}$", 
    "$C_{P_i}$", 
    ['Theory', 'BEMT'])
axs[1, 2].set_xlim([0, 0.01])

# Adjust layout to prevent overlapping
plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.show()


#######################################################################################################
# Problem 2
# Define the constants
sigma = 0.1
CT = 0.008
Nb2 = 2
Nb4 = 4
Cl_alpha = 2 * np.pi
theta_0 = (6 * CT / (sigma * Cl_alpha)) + 3/4 * sqrt(CT / 2)

# Define the range for r
r = np.linspace(0.01, 1, 50)

# Initialize lists for results
induced_lamda_2blade = []
induced_lamda_4blade = []
deltaCT_r_list_noloss = []
deltaCT_r_list_Nb2 = []
deltaCT_r_list_Nb4 = []
deltaCq_r_list_noloss = []
deltaCq_r_list_Nb2 = []
deltaCq_r_list_Nb4 = []
Cl_list_noloss = []
Cl_Nb2_list = []
Cl_Nb4_list = []

# Define the lambda function for sigma
sigma_function = lambda r: sigma

# Compute values for 2-blade and 4-blade configurations
lamda_no_loss_ideal = []
for r_i in r:
    # Initial guess for induced lambda
    induced_lamda_02 = induced_lamda_04 = sigma * Cl_alpha / (16 * 1) * (sqrt(1 + (32 * theta_0 * r_i / (sigma * Cl_alpha))) - 1)
    # Calculate no-tip-loss lambda for an ideal rotor
    lamda_ideal = sigma * Cl_alpha / 16 * (sqrt(1 + ((32 * theta_0 * r_i) / (sigma * Cl_alpha))) - 1)
    lamda_no_loss_ideal.append(lamda_ideal)

    # Iteratively compute induced lambda with tip loss factor
    for _ in range(4):
        F_2 = F(Nb2, r_i, induced_lamda_02)
        F_4 = F(Nb4, r_i, induced_lamda_04)
        induced_lamda_2 = sigma * Cl_alpha / (16 * F_2) * (sqrt(1 + (32 * F_2 * theta_0 * r_i / (sigma * Cl_alpha))) - 1)
        induced_lamda_4 = sigma * Cl_alpha / (16 * F_4) * (sqrt(1 + (32 * F_4 * theta_0 * r_i / (sigma * Cl_alpha))) - 1)
        induced_lamda_02 = induced_lamda_2
        induced_lamda_04 = induced_lamda_4

    # Store the final values for induced lambda
    induced_lamda_2blade.append(induced_lamda_02)
    induced_lamda_4blade.append(induced_lamda_04)

    # Compute other quantities without tip loss and for Nb2 and Nb4
    lamda_no_loss = sigma * Cl_alpha / 16 * (sqrt(1 + ((32 * theta_0 * r_i) / (sigma * Cl_alpha))) - 1)
    deltaCT_r_noloss = sigma * Cl_alpha / 2 * (theta_0 * r_i**2 - lamda_no_loss * r_i)
    deltaCT_r_Nb2 = sigma * Cl_alpha / 2 * (theta_0 * r_i**2 - induced_lamda_02 * r_i)
    deltaCT_r_Nb4 = sigma * Cl_alpha / 2 * (theta_0 * r_i**2 - induced_lamda_04 * r_i)
    
    deltaCq_r_noloss = 1/2 * sigma * r_i**2 * lamda_no_loss * Cl_alpha * (theta_0 - lamda_no_loss / r_i)
    deltaCq_r_Nb2 = 1/2 * sigma * r_i**2 * induced_lamda_02 * Cl_alpha * (theta_0 - induced_lamda_02 / r_i)
    deltaCq_r_Nb4 = 1/2 * sigma * r_i**2 * induced_lamda_04 * Cl_alpha * (theta_0 - induced_lamda_04 / r_i)
    
    Cl_noloss = Cl_alpha * (theta_0 - lamda_no_loss / r_i)
    Cl_Nb2 = Cl_alpha * (theta_0 - induced_lamda_02 / r_i)
    Cl_Nb4 = Cl_alpha * (theta_0 - induced_lamda_04 / r_i)
    
    deltaCT_r_list_noloss.append(deltaCT_r_noloss)
    deltaCT_r_list_Nb2.append(deltaCT_r_Nb2)
    deltaCT_r_list_Nb4.append(deltaCT_r_Nb4)
    
    deltaCq_r_list_noloss.append(deltaCq_r_noloss)
    deltaCq_r_list_Nb2.append(deltaCq_r_Nb2)
    deltaCq_r_list_Nb4.append(deltaCq_r_Nb4)
    
    Cl_list_noloss.append(Cl_noloss)
    Cl_Nb2_list.append(Cl_Nb2)
    Cl_Nb4_list.append(Cl_Nb4)

# Correct the last values to avoid division by zero
deltaCT_r_list_Nb2[-1] = 0
deltaCT_r_list_Nb4[-1] = 0
deltaCq_r_list_Nb2[-1] = 0
deltaCq_r_list_Nb4[-1] = 0
Cl_Nb2_list[-1] = 0
Cl_Nb4_list[-1] = 0

# Plotting the results
plt.figure(figsize=(12, 8))

# Plot for induced lambda
plt.subplot(2, 2, 1)
plt.plot(r, induced_lamda_2blade, 'x', linewidth=1.2, label='2 blades')
plt.plot(r, induced_lamda_4blade, 'x', linewidth=1.2, label='4 blades')
plt.plot(r, lamda_no_loss_ideal, '-', linewidth=1.2, label='no tip loss')
plt.title("Question 2: f(r) = λ")
plt.ylabel("λ")
plt.xlabel("r")
plt.legend()

# Plot for deltaCT_r
plt.subplot(2, 2, 2)
plt.plot(r, deltaCT_r_list_Nb2, 'x', linewidth=1.2, label='2 blades')
plt.plot(r, deltaCT_r_list_Nb4, 'x', linewidth=1.2, label='4 blades')
plt.plot(r, deltaCT_r_list_noloss, '-', linewidth=1.2, label='no tip loss')
plt.title("Question 2: f(r) = C_T")
plt.ylabel("$dC_{T}$/dr")
plt.xlabel("r")
plt.legend()

# Plot for deltaCq_r
plt.subplot(2, 2, 3)
plt.plot(r, deltaCq_r_list_Nb2, 'x', linewidth=1.2, label='2 blades')
plt.plot(r, deltaCq_r_list_Nb4, 'x', linewidth=1.2, label='4 blades')
plt.plot(r, deltaCq_r_list_noloss, '-', linewidth=1.2, label='no tip loss')
plt.title("Question 2: f(r) = C_Q")
plt.ylabel("$dC_{Q}$/dr")
plt.xlabel("r")
plt.legend()

# Plot for Cl
plt.subplot(2, 2, 4)
plt.plot(r, Cl_Nb2_list, 'x', label='2 blades')
plt.plot(r, Cl_Nb4_list, 'x', linewidth=1.2, label='4 blades')
plt.plot(r, Cl_list_noloss, '-', linewidth=1.2, label='no tip loss')
plt.title("Question 2: f(r) = C_L")
plt.ylabel("$C_l$")
plt.xlabel("r")
plt.legend()

plt.tight_layout()
plt.show()


######################################################################################
# Constants and parameters
Cl_alpha = 2 * np.pi
Cd0 = lambda alpha: 0.011
CT_req = 0.008
sigma = lambda r: 0.1
Nb = 4
N = 50
r_BEMT = np.linspace(1/N/2, 1-1/N/2, N)
theta_tw = np.arange(-30, 15, 5)  # degrees

# plt.figure(7)
# plt.suptitle('Problem 3: Effect of twist rate on radial distribution')

# for t_tw in theta_tw:
#     theta_0_n = theta_0_solver(True, CT_req, Nb, np.deg2rad(t_tw), sigma, Cl_alpha, Cd0)
#     theta_n = lambda r: theta_0_n + np.deg2rad(t_tw) * r
#     lambda_n, CT_n, delta_Cpi_t_n, Cp0_n, Cp_n, Cl_n, CTn, CPn = BEMT(True, N, Nb, sigma, theta_n, Cl_alpha, Cd0)

#     plt.subplot(2, 2, 1)
#     plt.plot(r_BEMT, lambda_n,'--', label=f'θ_tw={t_tw}°')
#     plt.ylabel("Induced inflow coefficient, λ")
#     plt.xlabel("r")
    
#     plt.subplot(2, 2, 2)
#     plt.plot(r_BEMT, CT_n,'--', label=f'θ_tw={t_tw}°')
#     plt.ylabel("Thrust distribution, $dC_T$/dr")
#     plt.xlabel("r")
    
#     plt.subplot(2, 2, 3)
#     plt.plot(r_BEMT, Cl_n,'--', label=f'θ_tw={t_tw}°')  
#     plt.ylabel("Lift coefficient, $C_L$")
#     plt.xlabel("r")
    
#     plt.subplot(2, 2, 4)
#     plt.plot(r_BEMT, delta_Cpi_t_n,'--', label=f'θ_tw={t_tw}°')
#     plt.ylabel("Inducced torque distribution, $dC_Q$/dr")
#     plt.xlabel("r")
 

plt.figure(figsize=(10, 6))
plt.suptitle('Problem 3: Effect of twist rate on radial distribution - Induced Inflow Coefficient')

for t_tw in theta_tw:
    theta_0_n = theta_0_solver(True, CT_req, Nb, np.deg2rad(t_tw), sigma, Cl_alpha, Cd0)
    theta_n = lambda r: theta_0_n + np.deg2rad(t_tw) * r
    lambda_n, CT_n, delta_Cpi_t_n, Cp0_n, Cp_n, Cl_n, CTn, CPn = BEMT(True, N, Nb, sigma, theta_n, Cl_alpha, Cd0)
    plt.plot(r_BEMT, lambda_n, '--', label=f'θ_tw={t_tw}°')

plt.ylabel("Induced inflow coefficient, $λ_i$")
plt.xlabel("r")
plt.legend()
plt.text(0.7, 0.9, r'$\sigma = 0.1, C_T = 0.008$', transform=plt.gca().transAxes) 
plt.show()


# Plot： Thrust distribution, dC_T/dr
plt.figure(figsize=(10, 6))
plt.suptitle('Problem 3: Effect of twist rate on radial distribution - Thrust Distribution')

for t_tw in theta_tw:
    theta_0_n = theta_0_solver(True, CT_req, Nb, np.deg2rad(t_tw), sigma, Cl_alpha, Cd0)
    theta_n = lambda r: theta_0_n + np.deg2rad(t_tw) * r
    lambda_n, CT_n, delta_Cpi_t_n, Cp0_n, Cp_n, Cl_n, CTn, CPn = BEMT(True, N, Nb, sigma, theta_n, Cl_alpha, Cd0)
    plt.plot(r_BEMT, CT_n, '--', label=f'θ_tw={t_tw}°')

plt.ylabel("Thrust distribution, $dC_T$/dr")
plt.xlabel("r")
plt.legend()
plt.text(0.7, 0.9, r'$\sigma = 0.1, C_T = 0.008$', transform=plt.gca().transAxes) 
plt.show()

# Plot： Lift coefficient, C_L
plt.figure(figsize=(10, 6))
plt.suptitle('Problem 3: Effect of twist rate on radial distribution - Lift Coefficient')

for t_tw in theta_tw:
    theta_0_n = theta_0_solver(True, CT_req, Nb, np.deg2rad(t_tw), sigma, Cl_alpha, Cd0)
    theta_n = lambda r: theta_0_n + np.deg2rad(t_tw) * r
    lambda_n, CT_n, delta_Cpi_t_n, Cp0_n, Cp_n, Cl_n, CTn, CPn = BEMT(True, N, Nb, sigma, theta_n, Cl_alpha, Cd0)
    plt.plot(r_BEMT, Cl_n, '--', label=f'θ_tw={t_tw}°')

plt.ylabel("Lift coefficient, $C_L$")
plt.xlabel("r")
plt.legend() 
plt.text(0.7, 0.9, r'$\sigma = 0.1, C_T = 0.008$', transform=plt.gca().transAxes, ha='center') 
plt.show()


# Plot： Induced torque distribution, dC_Q/dr
plt.figure(figsize=(10, 6))
plt.suptitle('Problem 3: Effect of twist rate on radial distribution - Induced Torque Distribution')

for t_tw in theta_tw:
    theta_0_n = theta_0_solver(True, CT_req, Nb, np.deg2rad(t_tw), sigma, Cl_alpha, Cd0)
    theta_n = lambda r: theta_0_n + np.deg2rad(t_tw) * r
    lambda_n, CT_n, delta_Cpi_t_n, Cp0_n, Cp_n, Cl_n, CTn, CPn = BEMT(True, N, Nb, sigma, theta_n, Cl_alpha, Cd0)
    plt.plot(r_BEMT, delta_Cpi_t_n, '--', label=f'θ_tw={t_tw}°')

plt.ylabel("Induced torque distribution, $dC_{pi}$/dr")
plt.xlabel("r")
plt.legend()
plt.text(0.7, 0.9, r'$\sigma = 0.1, C_T = 0.008$', transform=plt.gca().transAxes) 
plt.show()

   
# Plot: Effect of twist rate on rotor power versus thrust
plt.figure()
plt.suptitle('Problem 3: Effect of twist rate on rotor power versus thrust')
for t_tw in theta_tw:
    plt.plot(r_BEMT, Cp0_n,'--', label=f'θ_tw={t_tw}°')

plt.ylabel(" profile power distribution, $dC_{p0}/dr$")
plt.xlabel("r")
plt.text(0.7, 0.9, r'$\sigma = 0.1, N_b = 4$', transform=plt.gca().transAxes) 
# plt.xlim(0, 0.1)
# plt.ylim(0, 0.0001)
plt.show()


###############################################################################
# Problem 4
# Initial Parameters
Cl_alpha = 2 * np.pi
Cd_function = lambda alpha: 0.011 - 0.025 * alpha + 0.65 * alpha ** 2
sigma = lambda r: 0.1
Nb = 4

# Array Creation
N = 80
r_BEMT = np.linspace(0.01, 1, N)
# theta_tw = np.arange(-24, 1, 4)  # Adjusted to include 0
theta_tw = [-24, -20, -16, -12, -8, -4, 0]

# Plotting Setup
plt.figure(figsize=(10, 15))
plt.suptitle(f'Question 4: Power, Power factor (κ), FM\n#blades={Nb}')
lineStyles_q4 = ["*", "x", ">", "pentagram", "o", "square", "hexagram"]

# Creating subplots
ax1 = plt.subplot(3, 1, 1)
ax2 = plt.subplot(3, 1, 2)
ax3 = plt.subplot(3, 1, 3)


twist_angles = [-24, -20, -16, -12, -8, -4, 0]
labels = ['θ_tw=-24°', 'θ_tw=-20°', 'θ_tw=-16°', 'θ_tw=-12°', 'θ_tw=-8°', 'θ_tw=-4°', 'θ_tw=0°']
marker_styles = ["o", "*", ">", "s", "D", "x", "^"]

# Main Loop
for i, twist_ang in enumerate(theta_tw):
    CT_list = np.arange(0.002, 0.0105, 0.0005)
    CP_list = np.zeros(len(CT_list))
    Cpi_list = np.zeros(len(CT_list))
    
    for n in range(len(CT_list)):
        CT_req = CT_list[n]
        theta_0_n = theta_0_solver(True, CT_req, Nb, np.deg2rad(twist_ang), sigma, Cl_alpha, Cd_function)
        theta_n = lambda r: theta_0_n + np.deg2rad(twist_ang) * r
        _, _, delta_Cpi_t, _, _, _, _, CPn = BEMT(True, N, Nb, sigma, theta_n, Cl_alpha, Cd_function)
        CP_list[n] = CPn
        Cpi_list[n] = np.sum(delta_Cpi_t / N)
    
    k_factor = Cpi_list / (CT_list ** (3/2) / np.sqrt(2))
    FM = (CT_list ** (3/2) / np.sqrt(2)) / CP_list

    # Plotting with labels and legends, using different marker styles
    ax1.plot(CT_list, CP_list, linestyle="none", marker=marker_styles[i], label=f'θ_tw={twist_ang}°')
    ax2.plot(CT_list, k_factor, linestyle="none", marker=marker_styles[i], label=f'θ_tw={twist_ang}°')
    ax3.plot(CT_list, FM, linestyle="--", marker=marker_styles[i], label=f'θ_tw={twist_ang}°')

# Adjusting and labeling each subplot
ax1.set_title("f(C_P) = C_T")
ax1.set_ylabel("C_P (hover)")
ax1.set_xlabel("C_T")
ax1.set_xlim([0, 0.01])
ax1.legend()

ax2.set_title("f(κ) = C_T")
ax2.set_ylabel("κ (Induced Power Factor)")
ax2.set_xlabel("C_T")
ax2.set_xlim([0, 0.01])
ax2.set_ylim([1, 1.35])
ax2.legend()

ax3.set_title("f(FM) = C_T")
ax3.set_ylabel("FM")
ax3.set_xlabel("C_T")
ax3.set_xlim([0, 0.01])
ax3.legend()

plt.tight_layout()
plt.show()


print('Induced power factor kappa = {}'.format(k_factor))
###############################################################################
# Problem 5
# Constants
Cl_alpha = 2 * pi
Cd_alpha = lambda alpha: 0.008
Nb = 2
theta_tw = np.deg2rad(-15)  # Convert to radians
taper_ratio = np.array([1, 2, 3])  # σ_0:σ_tip = 1:1 to 3:1  
sigma_weighted = 0.1  # Thrust Weighted Solidity
CT_req = 0.008
N = 80
r_BEMT = np.linspace(0.01, 1, N)


# Plot generation
fig, axs = plt.subplots(2, 2, figsize=(12, 8))
plt.suptitle(f'Question 5: #blades = {Nb}; C_T={CT_req}')

for gama in taper_ratio:
    sigma_initial = sigma_weighted / (1 + 3/4 * (1/gama - 1))
    sigma_iter = sigma_initial * (1/gama - 1)
    sigma_function = lambda r: sigma_initial + sigma_iter * r
    theta_0 = theta_0_solver(True, CT_req, Nb, theta_tw, sigma_function, Cl_alpha, Cd_alpha)
    theta_function = lambda r: theta_0 + theta_tw * r

    lambda_induced, delta_CT_t_, delta_Cpi_t_, delta_Cp0_t, delta_Cp_t, Cl, CT, CP = \
        BEMT(True, N, Nb, sigma_function, theta_function, Cl_alpha, Cd_alpha)

    line_style = '-' if gama == 1 else '--'  # Solid line for 1:1, dashed for others
    label = f'Taper Ratio = {gama}:1 (Rectangular)' if gama == 1 else f'Taper Ratio = {gama}:1'
    
    axs[0, 0].plot(r_BEMT, lambda_induced, label=label, linestyle=line_style, linewidth=1.2)
    axs[0, 1].plot(r_BEMT, delta_CT_t_, label=label, linestyle=line_style, linewidth=1.2)
    axs[1, 0].plot(r_BEMT, delta_Cpi_t_, label=label, linestyle=line_style, linewidth=1.2)
    axs[1, 1].plot(r_BEMT, Cl, label=label, linestyle=line_style, linewidth=1.2)

for ax in axs.flat:
    ax.set_xlabel("r ")
    ax.legend()
axs[0, 0].set_title("Inflow Distribution")
axs[0, 0].set_ylabel("λ ")
axs[0, 1].set_title("Thrust Distribution")
axs[0, 1].set_ylabel("dC_T/dr ")
axs[1, 0].set_title("Induced Torque Distribution")
axs[1, 0].set_ylabel("dC_{Q_i}/dr = dC_{P_i}/dr")
axs[1, 1].set_title("Lift Coeff. Distribution")
axs[1, 1].set_ylabel("C_L ")

plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.show()


##########################################################################################
# Range of taper ratios and twist angles
taper_r = np.arange(0.5, 35.5, 0.5)
ang_tw = np.arange(-25, 6, 1)  # Degrees

# Initialize FM matrix
FM_fig = np.zeros((len(ang_tw), len(taper_r)))

# Iterative Calculation
for ang in range(len(ang_tw)):
    twist_ang = ang_tw[ang]
    for b in range(len(taper_r)):
        gama_ = taper_r[b]
        sigma_initial = sigma_weighted / (1 + 3/4 * (1/gama_ - 1))
        sigma_iter = sigma_initial * (1/gama_ - 1)
        sigma_function = lambda r: sigma_initial + sigma_iter * r
        theta_0 = theta_0_solver(True, CT_req, Nb, radians(twist_ang), sigma_function, Cl_alpha, lambda alpha: 0.008)
        theta_function = lambda r: theta_0 + radians(twist_ang) * r
        _, _, _, _, _, _, CT, CP = BEMT(True, N, Nb, sigma_function, theta_function, Cl_alpha, lambda alpha: 0.008)
        FM = (CT**(3/2) / sqrt(2)) / CP if CP != 0 else 0
        FM_fig[ang, b] = FM
        
# Find Maximum FM and Corresponding Parameters
FM_max = np.max(FM_fig)
t_tw_opt_index, TR_opt_index = np.unravel_index(np.argmax(FM_fig), FM_fig.shape)
twist_opt = ang_tw[t_tw_opt_index]
taper_opt = taper_r[TR_opt_index]

print('Problem 5:')
print("Maximum FM:")
print(FM_max)
print('\n Optimal angle of twist {} deg'.format(twist_opt))
print('\n Optimal taper ratio {} (Max FM)'.format(taper_opt))


# 2D Contour Plot with Specific Colorbar Range
plt.figure()
X, Y = np.meshgrid(ang_tw, taper_r)

# Define levels for contour plot
levels = np.arange(0.66, 0.82 + 0.005, 0.005) 

contour = plt.contourf(X, Y, FM_fig.T, levels=levels, cmap='viridis')
plt.colorbar(contour)

max_FM_value = FM_fig[t_tw_opt_index, TR_opt_index]
plt.scatter(twist_opt, taper_opt, color='red', marker='*', s=30)  # 减小标记大小

# 将文本注释进一步向左上方移动
offset_x = 2.5  # 水平偏移量
offset_y = 2.5  # 垂直偏移量
plt.text(twist_opt - offset_x, taper_opt + offset_y, f'({twist_opt}, {taper_opt}, {max_FM_value:.3f})', 
         color='black', ha='left', va='bottom')


# Rotate plot by 180 degrees
plt.gca().invert_xaxis()
plt.gca().invert_yaxis()

plt.title('Figure of Merit Distribution')
plt.xlabel('Twist Angle θ_tw (degrees)')
plt.ylabel('Taper Ratio')
plt.show()

# 3D Plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
X, Y = np.meshgrid(ang_tw, taper_r)
ax.plot_surface(X, Y, FM_fig.T, cmap='viridis', edgecolor='none')
ax.scatter(twist_opt, taper_opt, FM_max, color='red', marker='p', s=50)  # Highlight Max FM Point
ax.contour3D(X, Y, FM_fig.T, 30, colors='black', linestyles='dashed')  # Contour Lines
ax.set_xlabel('Twist Angle θ_tw (degrees)')
ax.set_ylabel('Taper Ratio')
ax.set_zlabel('Figure of Merit (FM)')
plt.show()
















































