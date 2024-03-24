# -*- coding: utf-8 -*-
"""
Created on Sun Nov 12 15:11:33 2023
Aero 614 hw 4
@author: yanyan

1. BONE (70 points)
The goal of this assignment is to further investigate the changes that age and environmental loading 
(1G, microgravity, etc.) have on geometry and bone mineral density of bone (modeled as a tubular structure). 
You are to perform simulations using a computer spreadsheet or Matlab. 
Credit will be given for discussion, professional work, the process articulated, and illustrations.

Assumption: Long bones can be simulated as hollow tubes, which are subjected in life mainly tobending stresses. 
Generally, the moment arms do not change with age although load magnitudes do change.

Starting conditions for young adult (assume 20 years old) femur shafts:
Female outer radius, ro = 1.45 cm Male outer radius, ro = 1.65 cm
Female inner radius, ri = 1.00 cm Male inner radius, ri = 1.20 cm

Assume that remodeling normally increases inner radius at the rate of 0.005 cm/year (based on measured data). 
The geometry of a hollow tube with effective mineral density of solid bone ρm = 1.05 g/cm3

Questions:
1.1.(30 points) 
A) Find the section modulus of young adult femoral shafts (for males and females). Assume that this is maintained through life (assuming constant skeletal loading).
                                                           
B) Derive an expression for ro as a function of changing ri. 
                                                           
C) Show how the section modulus can be maintained through age 80 (plot, for males and females, the evolution of the femur shaft inner and outer radius over time, from 20 years old to 80 years old). 
                                                                  
D) What happens to the BMD? (Plot BMD over time, from 20 years old to 80 years old, for males and females).

1.2.(20 point) We are interested in spaceflight where skeletal loading is significantly altered.
We would like to investigate the influence of higher remodeling rates. 

A) What happens to BMD (and other parameters) if remodeling rates (increase in ri) double? (Include updated versions of the previous plots in your response). 

B) If remodeling rates quadruple? (Include updated versions of the previous plots in your response). Assume constant section modulus (Is this a reasonable assumption?).

1.3.(20 point) Lets be more realistic and admit that the assumption of constant skeletal loading
is unrealistic. Investigate a loading profile that decreases linearly between ages 35 and 60,
being reduced by 30% at age at age 60, remaining constant thereafter (see loading profile
below). 

A) Explicitly explain/justify your loading profile(s) (Do you think this is a reasonable profile?). 
                                                       
B) Assume a linear effect on section modulus. Describe and show the results for BMD and other parameters? 
(Include updated version of plots for section modulus, femur shaft inner and outer radii, 
 and BMD over time from 20 years old to 80 years old, for males and females).

"""
import numpy as np
import matplotlib.pyplot as plt


def inertia(ro, ri):            # I = cross sectional moment of Inertia
     I = np.pi/4*(ro**4 - ri**4)
     return I


def section_modulus(ro, ri):
    z =  np.pi/4*(ro**4 - ri**4) / ro
    return z


# Initial conditions:
ro_female = 1.45    # outter
ro_male = 1.65
ri_female = 1       # inner
ri_male = 1.2

z_female = section_modulus(ro_female, ri_female)
z_male = section_modulus(ro_male, ri_male)
print('1.1a :')
print('Section modulus of female Z = {}'.format(z_female))
print('Section modulus of male Z = {}'.format(z_male))
print('\n1.1b:')
print('Z = (0.25*π*ro^4 - Z*ro = 0.25*π*(ri_initial + 0.005t)^4 ')


# Function to calculate the section modulus
def section_modulus(ro, ri):
    return 0.25 * np.pi * (ro**4 - ri**4) / ro

# Function to simulate the bone growth over time and plot the results
def simulate_bone_growth(ro_init, ri_init, ri_rate, rho, gender, years=70):
    # Calculate the initial section modulus
    Z = section_modulus(ro_init, ri_init)
    
    # Initialize arrays to store the results
    ages = np.arange(25, 25 + years)
    ro_values = np.full_like(ages, ro_init, dtype=float)
    ri_values = ri_init + ri_rate * (ages - 25)
    BMD_values = np.empty_like(ages, dtype=float)
    
    # Perform the calculations over the years
    for i, age in enumerate(ages):
        ro = ro_values[i]
        ro_side = 0.25 * np.pi * ro**4 - Z * ro  # one side of the equation for Z
        ri = ri_values[i]
        ri_side = 0.25 * np.pi * ri**4  # other side of the equation for Z
        
        # Iteratively increase ro until the equation is satisfied
        while ro_side <= ri_side:
            ro += 0.00005
            ro_side = 0.25 * np.pi * ro**4 - Z * ro
        
        # Update the outer radius for the current year
        ro_values[i] = ro
        
        # Calculate the area and BMD once the two sides of the equation converge
        area = np.pi * (ro**2 - ri**2)
        BMD_values[i] = area * rho / (2 * ro)
    
    # Plotting the results
    plt.figure(figsize=(10, 6))
    plt.plot(ages, ro_values, 'o-', label='Outer Radius')
    plt.plot(ages, ri_values, 's-', label='Inner Radius')
    plt.title('{}: Cortical Bone Radii Changes with Age'.format(gender))
    plt.xlabel('Age (years)')
    plt.ylabel('Radius (cm)')
    plt.legend()
    plt.grid(True)
    plt.show()
    
    return ages, ro_values, ri_values, BMD_values



rho = 1.05  # g/cm^3, density of solid bone
ri_rate = 4 * 0.005  # cm/year, increased remodeling rate


# print('Simulating for female:')
ages_female, ro_values_female, ri_values_female, BMD_values_female = simulate_bone_growth(
    ro_female, ri_female, ri_rate, rho, 'Female')

# print('Simulating for male:')
ages_male, ro_values_male, ri_values_male, BMD_values_male = simulate_bone_growth(
    ro_male, ri_male, ri_rate, rho, 'Male')

########################################################################################
# 1.1 D: BMD
# Initialize constants and conditions
ri_rate = 0.005  # cm/year, rate of increase of inner radius
rho = 1.05  # g/cm^3, density of solid bone
years = 70  # years to simulate

# Initialize results matrices
bone_results_male = np.zeros((years, 7))
bone_results_female = np.zeros((years, 7))

Z_male = section_modulus(ro_male, ri_male)  # section modulus for males
Z_female = section_modulus(ro_female, ri_female)  # section modulus for females

# Function to simulate the changes in BMD over time
def simulate_bone_changes(ro_init, ri_init, ri_rate, Z, rho, years):
    bone_results = np.zeros((years, 2))
    for t in range(years):
        age = t + 25  # starting from age 25
        ri = ri_init + ri_rate * t
        ro = ro_init  # start with initial outer radius for iteration
        # Iteratively find ro to maintain constant section modulus
        while 0.25 * np.pi * ro**4 - Z * ro <= 0.25 * np.pi * ri**4:
            ro += 0.00005
        # Calculate BMD
        BMD = np.pi * (ro**2 - ri**2) * rho / (2 * ro)
        # Store the results
        bone_results[t, 0] = age
        bone_results[t, 1] = (BMD - np.pi * (ro_init**2 - ri_init**2) * rho / (2 * ro_init)) / (np.pi * (ro_init**2 - ri_init**2) * rho / (2 * ro_init)) * 100
    return bone_results

# Simulate for males and females
results_male = simulate_bone_changes(ro_male, ri_male, ri_rate, Z_male, rho, years)
results_female = simulate_bone_changes(ro_female, ri_female, ri_rate, Z_female, rho, years)

# Plotting the results
plt.figure(figsize=(12, 8))
plt.plot(results_male[:, 0], results_male[:, 1], 's-', label='Males')
plt.plot(results_female[:, 0], results_female[:, 1], 'o-', label='Females')
plt.title('Variation in Bone Mineral Density Over Time with Consistent Skeletal Loading and a Remodeling Rate of 0.005 cm/year')
plt.xlabel('Age (years)')
plt.ylabel('% Change from starting conditions')
plt.legend()
plt.grid(True)
plt.show()


########################################################################################
# Function to calculate the section modulus for a given outer and inner radius
def section_modulus(ro, ri):
    return 0.25 * np.pi * (ro**4 - ri**4) / ro

def simulate_bone_changes(ro_init, ri_init, ri_rate, rho, years=70):
    Z = section_modulus(ro_init, ri_init)  # Calculate initial section modulus
    age = np.arange(25, 25 + years)  # Age starts from 25
    ro_change = np.zeros(years)  # Percentage change in outer radius
    ri_change = np.zeros(years)  # Percentage change in inner radius
    BMD_change = np.zeros(years)  # Percentage change in BMD

    # Initial BMD needs to be calculated inside the function
    BMD_initial = np.pi * (ro_init**2 - ri_init**2) * rho / (2 * ro_init)

    for t in range(years):
        ri = ri_init + ri_rate * t  # Inner radius increases each year
        ro = ro_init  # Outer radius is recalculated each year to maintain Z

        # Iteratively find the new outer radius to maintain section modulus
        while True:
            ro_side = 0.25 * np.pi * ro**4 - Z * ro  # Calculate one side of the equation for Z
            ri_side = 0.25 * np.pi * ri**4  # Calculate the other side of the equation for Z
            if ro_side < ri_side:
                ro += 0.00005  # Increment outer radius
            else:
                break
        
        # Recalculate BMD for the current year
        BMD = np.pi * (ro**2 - ri**2) * rho / (2 * ro)
        ro_change[t] = (ro - ro_init) / ro_init * 100  # Store % change in outer radius
        ri_change[t] = (ri - ri_init) / ri_init * 100  # Store % change in inner radius
        BMD_change[t] = (BMD - BMD_initial) / BMD_initial * 100  # Store % change in BMD

    return age, ro_change, ri_change, BMD_change



# Initial conditions for female and male
ro_init_f = 1.45  # cm, initial outer radius for females
ri_init_f = 1.0  # cm, initial inner radius for females
ro_init_m = 1.65  # cm, initial outer radius for males
ri_init_m = 1.20  # cm, initial inner radius for males


# Remodeling rates
doubled_rate = 0.005*2  # cm/year
quadrupled_rate = 0.005*4  # cm/year

# Run simulations
age, ro_change_f_doubled, ri_change_f_doubled, BMD_change_f_doubled = simulate_bone_changes(ro_init_f, ri_init_f, doubled_rate, rho)
age, ro_change_m_doubled, ri_change_m_doubled, BMD_change_m_doubled = simulate_bone_changes(ro_init_m, ri_init_m, doubled_rate, rho)
age, ro_change_f_quadrupled, ri_change_f_quadrupled, BMD_change_f_quadrupled = simulate_bone_changes(ro_init_f, ri_init_f, quadrupled_rate, rho)
age, ro_change_m_quadrupled, ri_change_m_quadrupled, BMD_change_m_quadrupled = simulate_bone_changes(ro_init_m, ri_init_m, quadrupled_rate, rho)

# Create figure and axes for the 2x2 grid
fig, axs = plt.subplots(2, 2, figsize=(14, 10), constrained_layout=True)

# Plot results for females with doubled remodeling rate
axs[0, 0].plot(age, ri_change_f_doubled, 'o-', label='Inner Radius Change')
axs[0, 0].plot(age, ro_change_f_doubled, 's-', label='Outer Radius Change')
axs[0, 0].plot(age, BMD_change_f_doubled, '^-', label='BMD Change')
axs[0, 0].set_title('Female with 0.01 cm/year remodeling')
axs[0, 0].set_xlabel('Age (years)')
axs[0, 0].set_ylabel('% Change')
axs[0, 0].legend()
axs[0, 0].grid(True)

# Plot results for males with doubled remodeling rate
axs[0, 1].plot(age, ri_change_m_doubled, 'o-', label='Inner Radius Change')
axs[0, 1].plot(age, ro_change_m_doubled, 's-', label='Outer Radius Change')
axs[0, 1].plot(age, BMD_change_m_doubled, '^-', label='BMD Change')
axs[0, 1].set_title('Male with 0.01 cm/year remodeling')
axs[0, 1].set_xlabel('Age (years)')
axs[0, 1].set_ylabel('% Change')
axs[0, 1].legend()
axs[0, 1].grid(True)

# Plot results for females with quadrupled remodeling rate
axs[1, 0].plot(age, ri_change_f_quadrupled, 'o-', label='Inner Radius Change')
axs[1, 0].plot(age, ro_change_f_quadrupled, 's-', label='Outer Radius Change')
axs[1, 0].plot(age, BMD_change_f_quadrupled, '^-', label='BMD Change')
axs[1, 0].set_title('Female with 0.02 cm/year remodeling')
axs[1, 0].set_xlabel('Age (years)')
axs[1, 0].set_ylabel('% Change')
axs[1, 0].legend()
axs[1, 0].grid(True)

# Plot results for males with quadrupled remodeling rate
axs[1, 1].plot(age, ri_change_m_quadrupled, 'o-', label='Inner Radius Change')
axs[1, 1].plot(age, ro_change_m_quadrupled, 's-', label='Outer Radius Change')
axs[1, 1].plot(age, BMD_change_m_quadrupled, '^-', label='BMD Change')
axs[1, 1].set_title('Male with 0.02 cm/year remodeling')
axs[1, 1].set_xlabel('Age (years)')
axs[1, 1].set_ylabel('% Change')
axs[1, 1].legend()
axs[1, 1].grid(True)

plt.show()

#################################################################################


# Function to define the loading profile
def loading_profile(age):
    if age < 35:
        return 1.0  # 100% loading
    elif 35 <= age <= 60:
        # Linear decrease from 100% to 70% over 25 years
        return 1 - 0.012 * (age - 35)
    else:
        return 0.7  # 70% loading after age 60

# Function to simulate bone changes with variable skeletal loading
def simulate_bone_changes_variable_loading(ro_init, ri_init, ri_rate, rho, years=70, Z_init=None):
    if Z_init is None:
        Z_init = 0.25 * np.pi * (ro_init**4 - ri_init**4) / ro_init
    age = np.arange(25, 25 + years)
    ro_change = np.zeros(years)
    ri_change = np.zeros(years)
    BMD_change = np.zeros(years)
    BMD_initial = np.pi * (ro_init**2 - ri_init**2) * rho / (2 * ro_init)
    
    for t in range(years):
        current_age = age[t]
        ri = ri_init + ri_rate * t
        ro = ro_init
        Z = Z_init * loading_profile(current_age)
        
        # Find the outer radius that corresponds to the current loading profile
        while 0.25 * np.pi * ro**4 - Z * ro <= 0.25 * np.pi * ri**4:
            ro += 0.00005
        
        # Calculate the area and BMD for the current year
        area = np.pi * (ro**2 - ri**2)
        BMD = area * rho / (2 * ro)
        
        # Store the results as percentage changes from the initial values
        ro_change[t] = ((ro - ro_init) / ro_init) * 100
        ri_change[t] = ((ri - ri_init) / ri_init) * 100
        BMD_change[t] = ((BMD - BMD_initial) / BMD_initial) * 100
    
    return age, ro_change, ri_change, BMD_change




# Remodeling rate
remodeling_rate = 0.005  # cm/year

# Run simulations for variable loading
age, ro_change_f_variable, ri_change_f_variable, BMD_change_f_variable = simulate_bone_changes_variable_loading(ro_init_f, ri_init_f, remodeling_rate, rho)
age, ro_change_m_variable, ri_change_m_variable, BMD_change_m_variable = simulate_bone_changes_variable_loading(ro_init_m, ri_init_m, remodeling_rate, rho)

# Create figure and axes for the plots
fig, axs = plt.subplots(2, figsize=(10, 12))

# Plot results for female with variable skeletal loading
axs[0].plot(age, ri_change_f_variable, 'o-', label='Inner Radius Change')
axs[0].plot(age, ro_change_f_variable, '--', label='Outer Radius Change')
axs[0].plot(age, BMD_change_f_variable, '^-', label='BMD Change')
axs[0].set_title('Female with variable skeletal loading and remodeling rate of 0.005 cm/year')
axs[0].set_xlabel('Age (years)')
axs[0].set_ylabel('% Change from baseline conditions')
axs[0].legend()
axs[0].grid(True)

# Plot results for male with variable skeletal loading
axs[1].plot(age, ri_change_m_variable, 'o-', label='Inner Radius Change')
axs[1].plot(age, ro_change_m_variable, '--', label='Outer Radius Change')
axs[1].plot(age, BMD_change_m_variable, '^-', label='BMD Change')
axs[1].set_title('Male with variable skeletal loading and remodeling rate of 0.005 cm/year')
axs[1].set_xlabel('Age (years)')
axs[1].set_ylabel('% Change from baseline conditions')
axs[1].legend()
axs[1].grid(True)

plt.show()







