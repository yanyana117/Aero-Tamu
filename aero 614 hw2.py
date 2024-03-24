# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 21:45:26 2023

@author: yanyan
"""
import numpy as np
import matplotlib. pyplot as plt
import math
from scipy.optimize import curve_fit
import sympy as sp

'''
Calculate ESPVR line:
    point 1 (15, 0)
    point 2 (70,102)
    point 3 (40,60)
    
    

'''
xx = np.array([ 15, 70,40 ])

yy = np.array([ 0, 102, 60])

# Find induced power factor (κ):lease square fit to measured data 
# 定义线性模型函数
def linear_model(x, m, b):
    return m * x + b

# 使用最小二乘法拟合线性模型
params, covariance = curve_fit(linear_model, xx, yy)

# 从拟合结果中获取斜率和截距
m_fit, b_fit = params

x = xx

y = m_fit*x + b_fit
 
# 打印估计的斜率和截距
print(f"估计的斜率 m: {m_fit}")
print(f"估计的截距 b: {b_fit}")
print(f'y = {m_fit:.5f}x +{b_fit:.5f}')
print()




print(1.83956*xx +-22.64835)
