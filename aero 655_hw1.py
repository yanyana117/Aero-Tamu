
# -*- coding: utf-8 -*-

"""
Created on Wed Sep 13 01:54:05 2023

@author: yanyan

Problem 2:
    Measurements have been made of rotor performance at a fixed rotor speed for a series of blade pitch angles.

    The values of C_T that were measured were 6.0000E-06, 0.0010490,0.0023760, 0.0040760 and 0.0055810, 
 and the corresponding values of CP were 0.000197, 0.000226, 0.000282, 0.000405 and 0.000555, respectively.

     Plot this data in the form of a power polar (CT vs. CP). Explain (and show in a chart) 
 how to extract induced power factor (κ) and zero thrust power (profile power) for the rotor from these measured data.
     
     Then, to the experimental power polar chart add the analytical power polar curve predicted by modified momentum theory. 
 Note: Plot all the charts in Matlab or Excel. Denote experimental data by points (no line) and analytical results by continuous line (no markers)]

    Translate:
    已经对一系列叶片桨距角进行了固定转子速度下的性能测量经对一系列叶片桨距角进行了固定转子速度下的性能测量。
    请以功率极坐标图的形式绘制这些数据（CT vs. CP）。
    解释并在图表中展示如何从这些测量数据中提取叶片桨的 感应功率因子（κ） 和 零推力功率（轮廓功率）。
    然后，在实验功率极坐标图上添加由改进的动量理论预测的分析功率极坐标曲线。
    [注意：在Matlab或Excel中绘制所有图表。用点表示实验数据（无线条），用连续线表示分析结果（无标记）]。

"""


import numpy as np
import matplotlib. pyplot as plt
import math
from scipy.optimize import curve_fit
import sympy as sp
# Given:
CT = np.array([ 6.0000E-06, 0.0010490, 0.0023760, 0.0040760, 0.0055810])
Cp = np.array([ 0.000197, 0.000226, 0.000282, 0.000405 , 0.000555])

# Simple momentum theory:
Cp_ideal = CT**(3/2)/math.sqrt(2)

# Set print options to display in scientific notation
np.set_printoptions(formatter={'float': '{:e}'.format})

# print('Cp_ideal = {}'.format(Cp_ideal))

# Calculate zero thrust power (profile power)
# print('Cp[0] = {}'.format(Cp[0]))


'''
xx = Cp_ideal
yy = Cp

y_sum2 = sum([yi**2 for yi in yy])
x_sum2 = sum([xi**2 for xi in xx])

xy_sum=sum([xi*yi for xi,yi in zip(xx,yy)])

x_sum = sum(xx)
y_sum = sum(yy)

print('sum(y^2) = {}'.format(y_sum2))
print('sum(x^2) = {}'.format(x_sum2))
print('sum(xy) = {}'.format(xy_sum))
print('sum(x) = {}'.format(x_sum))
print('sum(y) = {}'.format(y_sum))


from sympy import symbols, Eq, solve

# 定义符号变量
m, c = symbols('m c')

# 创建方程组
eq1 = Eq(2*x_sum2*m + 10*x_sum*c, 10*xy_sum)
eq2 = Eq(10*c + 10 * x_sum*m, 10*y_sum)

# 解方程组
solution = solve((eq1, eq2), (m, c))
print(solution)
'''




##################################################################
# Find induced power factor (κ):lease square fit to measured data 
# 定义线性模型函数
def linear_model(x, m, b):
    return m * x + b

# 使用最小二乘法拟合线性模型
params, covariance = curve_fit(linear_model, Cp_ideal, Cp)

# 从拟合结果中获取斜率和截距
m_fit, b_fit = params

x = Cp_ideal

y = m_fit*x + b_fit
 
# 打印估计的斜率和截距
print(f"估计的斜率 m: {m_fit}")
print(f"估计的截距 b: {b_fit}")
print(f'y = {m_fit:.5f}x +{b_fit:.5f}')
print()
##################################################################
# Check:
# y = mx +b
# Cp = m * Cp_ideal + Cp_0

Cp_imean = np.mean(Cp_ideal)
Cp_mean = np.mean(Cp)

func_up = np.sum( (Cp_ideal - Cp_imean)*(Cp - Cp_mean))
func_below = np.sum((Cp_ideal-Cp_imean)**2)    

m = func_up/func_below
Cp_0 = np.mean(Cp) - m*np.mean(Cp_ideal)


print('kappa = {}'.format(m))
print('Cp_0= {}'.format(Cp_0))


# Modified theory:
Cp_mt =  CT**(3/2)/math.sqrt(2) + Cp_0

# plot 1
plt.plot(CT,Cp_ideal, linestyle=':', label = 'Simple momentum theory')
plt.plot( CT, y, label='Modified momentum theory (κ)',color='red')
plt.plot( CT, Cp_mt, linestyle='--', label='Modified  theory')

plt.scatter(CT, Cp,label  = 'Measured')   # measured experimental data by points


plt.xlabel('Thrust coefficient, C_T')
plt.ylabel('Power coefficient, Cp')
plt.legend()
plt.grid(True)

plt.show()

# plot 2
plt.plot(Cp_ideal, CT, linestyle=':', label='Simple momentum theory')
plt.plot(y, CT, label='Modified momentum theory (κ)',color='red')
plt.plot(Cp_mt, CT, linestyle='--', label='Modified theory')
plt.scatter(Cp, CT, label='Measured')  # measured experimental data by points

plt.xlabel('Power coefficient, Cp')
plt.ylabel('Thrust coefficient, C_T')
plt.legend()
plt.grid(True)

plt.show()


# plot 3
plt.plot(x, y, label='y = mx + b',color='red')

plt.scatter(Cp_ideal, Cp, label='Measured')  # measured experimental data by points

plt.xlabel('Cp measurement ')
plt.ylabel('Cp ideal')
plt.legend()
plt.grid(True)




plt.show()

##############################################################################
print()
print('probelm 3')

'''
Problem 3. A helicopter with a gross weight of 3,000 lb (1,360.5 kg), a main rotor radius of 13.2 ft
(4.0 m), a rotor tip speed of 680 ft/s (207.3 m/s), and has 275 hp (205 kW) delivered to the
main rotor shaft. For hovering at sea level conditions, compute: (a) the rotor disk loading,
(b) the ideal power loading, (c) the thrust and torque coefficients, (d) the figure of merit
and actual power loading. [Note: Recommend working this problem in SI units]

Problem 4. For the helicopter described in the previous question, the tail rotor radius is 2.3 ft (0.701
m) and the tail rotor is located 15.3 ft (4.66 m) from the main rotor shaft. Calculate the
thrust and power required by the tail rotor when helicopter is hovering at sea level. Assume
that the figure of merit of the tail rotor is 0.7.

'''
# Problem 3
# initnal condition:
mass = 1360.5   # kg
R = 4   # m
v_tip = 207.3   # m/s
power_actual = 205000 # W
rho = 1.225 # kg/m^3 @sea level


# Calculate:
weight = mass*9.81
T = weight
A = math.pi*R**2 
DL = weight/A
power_ideal = T*np.sqrt(T/(2*rho*A))
PL_ideal = weight/power_ideal
Cq = power_actual/(rho*A*v_tip**3) 
Ct = T/(rho*(v_tip**2)*A)
FM = power_ideal/power_actual
PL_act = T/power_actual


# Problem 4:
# initnal condition:
R_tail = 0.701   # m
L_t2m = 4.66  # m; distance from tail to main rotor
FM_tail = 0.7


# Calculates:
w = v_tip/R    # angualr velocity
Q = power_actual/w  # torque
T_tail = Q/L_t2m
A_tail = math.pi*R_tail**2
v_tail = np.sqrt(T_tail/(2*rho*A_tail))
power_tail_ideal = T_tail*v_tail
power_tail_act = power_tail_ideal/FM_tail


print('weight = {} N '.format(T))
print('A = {} m^2 '.format(A))
print('DL = {} N/m^2 '.format(DL))
print('power_ideal = {} W '.format(power_ideal))
print('PL_ideal = {} N/W '.format(PL_ideal))
print('Cq torque coefficients = {:e} '.format(Cq))
print('Ct torque coefficients = {:e} '.format(Ct))
print('FM = {} '.format(FM))
print('PL_act = {} N/m '.format(PL_act))
print()
print('probelm 4')
print('Torque Q = {} Nm '.format(Q))
print('Thrust on tail = {} N '.format(T_tail))
print('Ideal Power on tail = {} W '.format(power_tail_ideal))
print('Actual Power on tail = {} W '.format(power_tail_act))




