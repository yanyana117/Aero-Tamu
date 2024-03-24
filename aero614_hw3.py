# -*- coding: utf-8 -*-
"""
Created on Mon Oct 23 19:34:15 2023

@author: yanyan
"""

import scipy.stats as stats
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm


def probability_between(mu, sigma, lower, upper):
    z_lower = (lower - mu) / sigma
    z_upper = (upper - mu) / sigma
    p_lower = norm.cdf(z_lower)
    p_upper = norm.cdf(z_upper)
    return p_upper - p_lower, z_lower,z_upper

def probability_less_than(mu, sigma, bound):
    z_bound = (bound - mu) / sigma
    return norm.cdf(z_bound)

def probability_larger_than(mu, sigma, bound):
    z_bound = (bound - mu) / sigma
    return (1 - norm.cdf(z_bound))

########################################################################################

# Define the range of values for the x-axis
x = np.linspace(-3, 3, 1000)


# Calculate the probability density function (PDF) for the standard normal distribution
pdf = stats.norm.pdf(x)


# Define the scenarios you want to calculate
scenarios = [
('a: P(0 < z < 2.07)', 0, 2.07, 'green'),
('b: P(-1.83 < z < 0)', -1.83, 0, 'blue'),
('c: P(-1.59 < z < 1.88)', -1.59, 1.88, 'orange'),
('d: P(1.33 < z < 1.88)', 1.33, 1.88, 'red'),
('e: P(-2.56 < z < 0.37)', -2.56, 0.37, 'purple'),
('f: P(z < 1.66)', -3, 1.66, 'pink'),
('g: P(z < -2.03)', -3, -2.03, 'cyan'),
('h: P(z > -1.13)', -1.13, 3, 'magenta'),
('i: P(z < 1.93)', -3, 1.93, 'brown'),
('j: P(z > -1.77)', -1.77, 3, 'lime')
]


# Create subplots for each scenario
fig, axes = plt.subplots(5, 2, figsize=(10, 10))


for (label, a, b, color), ax in zip(scenarios, axes.ravel()):
    ax.plot(x, pdf, label='Standard Normal Distribution')

    # Determine the type of probability we're dealing with
    if a is None:
        a = -np.inf  # 将None替换为-np.inf
        prob = probability_less_than(0, 1, b)
        fill_between_condition = x <= b
    elif b is None:
        b = np.inf  # 将None替换为np.inf
        prob = probability_larger_than(0, 1, a)
        fill_between_condition = x >= a
    else:
        prob,_,_ = probability_between(0, 1, a, b)
        fill_between_condition = (x >= a) & (x <= b)

    ax.fill_between(x, pdf, where=fill_between_condition, alpha=0.5, label=f'{label} = {prob:.4f}', color=color)
    ax.set_xlabel('Z')
    ax.set_ylabel('Probability Density')
    ax.set_title(label)
    ax.legend()
    ax.grid(True)

# Adjust spacing between subplots
plt.tight_layout()

# Show the figure with all the subplots
plt.show()

###########################################################################################

p2a,_,_  = probability_between(mu=27989, sigma=3250, lower=20000, upper=30000)
p2b = probability_less_than(mu=27989, sigma=3250, bound=20000)


p4a,z4a_low,z4a_up  = probability_between(mu=36.2, sigma=3.7, lower=36, upper=37.5)


print(f"2a: P = {p2a:.5f}")
print(f"2b: P = {p2b:.5f}")
print(f"4a: P = {p4a:.5f}")
print(f"4a: z_low = {z4a_low:.3f}")
print(f"4a: z_up = {z4a_up:.3f}")


###########################################################
# Define the range of values for the x-axis
x = np.linspace(-3, 3, 1000)

# Calculate the probability density function (PDF) for the standard normal distribution
pdf = stats.norm.pdf(x)

# Define the scenarios you want to calculate
scenarios = [
    ('2a: P(-2.46 < z < 0.62)', -2.46, 0.62, 'green'),
    ('2b: P(z < -2.46)', None, -2.46, 'red'),  # None表示负无穷
    ('3: P(z > 1.04)', 1.04, None, 'orange'),  # 同样，正无穷也是None
    (f'4a: P({z4a_low:.4f} < z < {z4a_up:.4f})', -0.05, 0.35, 'orange'),
    ('4b: P(-0.21 < z < 1.36)', -0.21, 1.36, 'orange'),
    # ('5: P(-0.21 < z < 1.36)', 1.96, None, 'orange')
    
]


# Create subplots for each scenario with 2 rows and 5 columns layout
fig, axes = plt.subplots(5, 2, figsize=(10,10))  # 5列2行

# Flatten the axes to make it easy to iterate over
axes = axes.flatten()

for (label, a, b, color), ax in zip(scenarios, axes):
    ax.plot(x, pdf, label='Standard Normal Distribution')
    if a is None:
        a = -np.inf  # 将None替换为-np.inf
        prob = stats.norm.cdf(b)
        fill_between_condition = x <= b
    elif b is None:
        b = np.inf  # 将None替换为np.inf
        prob = 1 - stats.norm.cdf(a)
        fill_between_condition = x >= a
    else:
        prob = stats.norm.cdf(b) - stats.norm.cdf(a)
        fill_between_condition = (x >= a) & (x <= b)
    
    ax.fill_between(x, pdf, where=fill_between_condition, alpha=0.5, label=f'{label} = {prob:.4f}', color=color)
    ax.set_xlabel('Z')
    ax.set_ylabel('Probability Density')
    ax.set_title(label)
    ax.legend()
    ax.grid(True)

# Adjust spacing between subplots
plt.tight_layout()

# Show the figure with all the subplots
plt.show()
