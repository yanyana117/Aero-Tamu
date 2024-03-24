# -*- coding: utf-8 -*-
"""
Created on Mon Oct 23 00:32:42 2023

@author: yanyan
"""


import numpy as np
import matplotlib. pyplot as plt
import math


# Given:
W = 6000    # lbf
vc =  np.array([-600, -300,0, 300, 600 ])/60 # ft/min, climb velocity, unit convert to ft/s
vc =  np.array([600,-600,300,-300,0 ])/60 
vc = np.sort(vc)
R = 21  # ft
FM = 0.75

rho_sea = 2.37717E-3   # slug/ft3
rho_5k = 2.04834E-3  # slug/ft3
rho_10k = 1.75549E-3 # slug/ft3
rho = np.array([rho_sea, rho_5k,rho_10k])
alitude = np.array(["Sea Level", "5000 ft", "10000 ft"])

# Calculate:
Thrust = W
A = math.pi*R**2
# A_book = math.pi*20**2

# power ratio
# if vc/vh > 0 :
#  P =
k = 1
k1= -1.125
k2= -1.372
k3= -1.78
k4= -0.655


# Initialize Output Lists: 初始化输出列表
"""
 在整个代码中，这些变量（output_rho、output_vcs 和 output_power）用于存储不同条件下的计算结果
"""

output_rho = [] # 用于存储不同海拔条件下的空气密度,在后续部分，循环遍历不同的密度值（对应不同的海拔条件），并将这些密度值添加到output_rho 列表中
output_vcs = [] # 存储不同爬升速度；   eg,这是就好像你准备了一个空的购物袋（output_vcs 列表）
output_power = [] # 这是用于存储不同情况下计算的功率值
output_vi = []
output_vivh = []
output_vcvh = []
output_vh=[]

# Calculation:
    
for density in rho: 
    vh = np.sqrt(Thrust / (2*density*A))    
  
                            
    for vc_current in vc:                                        
        output_vcs.append(vc_current)  
                                        
                                                                         
        if vc_current / vh >= 0:             
            """Use vi/vh equation for positive vc"""
            """Power = T*(vc+vi)"""
            vi = -vc_current/2 + np.sqrt((vc_current/2)**2 + vh**2)
            P_act = Thrust*(vc_current + vi)/0.75
            vivh = vi/vh
            vcvh = vc_current / vh
            
        elif vc_current / vh == 0:
            """Power = T*vh"""
            P_act = Thrust* vh/0.75
            
            
        elif -2 < vc_current / vh < 0:
            """Use k1 k2 k3 k4 equation"""
            vcvh = vc_current / vh
            vi = vh * (k+ k1*vcvh + k2 * vcvh**2 + k3 * vcvh**3 + k4 * vcvh**4)
            vivh = vi/vh
            
            """Power = T*(vc+vi)"""
            P_act = Thrust*(vc_current + vi)/0.75

        elif vc_current / vh <= -2:
            """Use vi/vh equation for negative vc"""
            vi = -vc_current/2 - np.sqrt(vc_current**2/4 - vh**2)
            vivh = vi/vh
            vcvh = vc_current / vh
            """Power = T*(vc+vi)"""
            P_act = Thrust*(vc_current + vi)/0.75
            
            
        else:
            raise RuntimeError("You messed up dude")
            """You did something wrong if this ever executes"""
            
        # Add final power to output_power
        output_power.append(P_act)  
        output_vi.append(vi)  
        output_vh.append(vh)

        output_vivh.append(vivh) 
        output_vcvh.append(vcvh)


print("y: vi/vh:")
print(output_vivh)
print('\n')

print("vivh_test:")
vivh_test = [vi/vh for vi, vh in zip(output_vi, output_vh) if vh != 0]
print(vivh_test)

print("x: vc/vh:")
print(output_vcvh)
print('\n')






# 使用matplotlib绘制图像
# xx = np.linspace(-4, 4, 15)
# plt.plot(output_vcs/vh , output_vi/vh, 'r-', label='Induced Velocity Ratio')
plt.plot(output_vcvh , output_vivh, 'b*', label='Induced Velocity Ratio')
plt.plot(np.array(output_vcs)/vh , np.array(output_vi)/vh, 'r-', label='Induced Velocity Ratio (Red Line)')
plt.grid(True)
plt.xlabel('Climb velocity ratio, Vc/Vh')
plt.ylabel('Induced velocity ratio, Vi/Vh')
plt.title('Induced velocity Curves')
plt.legend()
plt.show()













