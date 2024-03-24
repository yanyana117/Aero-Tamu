# -*- coding: utf-8 -*-
"""
Created on Mon Oct 16 23:24:54 2023
@author: yanyan
AERO 455/655 Helicopter Aerodynamics hw 2 (Fall 2023)

Problem 1:
Given a helicopter of weight, W = 6000 lb, calculate the power required in hover  and at rates of climb (Vc) of 300, 600, -300, and -600 ft/min.
The radius of main  rotor is 21 ft and the rotor has a figure of merit of 0.75.
Calculate your results at mean sea level conditions, and for density altitudes of five thousand and ten thousand feet.
Plot your results in the form of power  required versus climb velocity.
Discuss the factors that will determine the maximum vertical climb rate of a helicopter.

Probkem 2:
Starting from the momentum theory result for the induced flow at the rotor disk in forward flight, program numerical solutions to this inflow equation using both  the fixed-point (FP) iterative method and the Newton-Raphson (NR) method. 
Compare your numerical results with the exact analytic result for the case where rotor disk angle of attack (α) is zero. 
Plot analytic theory as a curve and numerical solution as discrete points with symbols. For both numerical methods, explore the effect of the initial guess for λ on the results from the numerical solution. [Assume CT = 0.01]

Problem 3:
Using the numerical solution in Q.2, plot some example graphs of inflow ratio and induced inflow ratio versus advance ratio in level flight for a series of rotor angles of attack (α). 
Use both positive and negative angles of attack (vary α from -4° to 8° in 2° increments). 
Explore the effect of initial guess for λ on the numerical solution, and suggest why a single unique initial condition may not always be adequate (use both FP and NR methods). 
Comment on your results.

Problem 4:
Compare the results from the iterative solution of the general inflow equation to the exact analytical inflow equation for the case of axial climbing and descending flight.
Numerically, reconstruct the complete induced velocity curve and the universal power curve, and compare with the exact (analytic) results derived in class. 
Plot the analytic theory as a curve and numerical solution as discrete points with symbols. Use both FP and NR methods to solve the inflow equation. 
Determine the conditions (if any) under which these methods may fail or give unexpected results.
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
k3= -1.718
k4= -0.655


# Initialize Output Lists: 初始化输出列表
"""
 在整个代码中，这些变量（output_rho、output_vcs 和 output_power）用于存储不同条件下的计算结果
"""

output_rho = [] # 用于存储不同海拔条件下的空气密度,在后续部分，循环遍历不同的密度值（对应不同的海拔条件），并将这些密度值添加到output_rho 列表中
output_vcs = [] # 存储不同爬升速度；   eg,这是就好像你准备了一个空的购物袋（output_vcs 列表）
output_power = [] # 这是用于存储不同情况下计算的功率值
output_vis = []
output_vivh = []
output_vcvh = []
output_vh=[]
# Calculation:
    
for density in rho: 
    vh = np.sqrt(Thrust / (2*density*A))    
    # density(current rho) 是在每次循环迭代中使用的临时变量，用于表示当前海拔条件下的空气密度。
    # rho (sum rho的合集)是一个包含了不同海拔条件下的多个密度值的列表，用于在循环中遍历不同的海拔条件。
    # 这是第一个循环：外部循环，它针对不同的密度值（对应不同的海拔条件）进行迭代的部分。
    # 在这个循环内，density(current rho) 是当前处理的密度值，它会依次取 rho 列表中的每个密度值，然后进行相关的计算。
    # 这个外部循环的目的是针对不同海拔条件下的密度值执行一组操作，然后将结果存储在 output_rho、output_vcs 和 output_power 列表中
                            
    for vc_current in vc:   # vc 是sum vc合集； 第二个内部循环，vc_current 是当前处理的爬升速度值，它会依次取 vc 数组中的每个速度值，并执行与不同爬升速度相关的操作。
                            # 内部循环是嵌套在外部循环内的，用于在不同的密度值条件下对各种爬升速度进行计算。
                            # eg, 我站在vc水果店前，准备挑选各种不同的水果（vc_current 代表不同的水果）；for vc_current in vc: 
                                       

        output_vcs.append(vc_current)   # output_vcs = []，output_vcs.是储存列表，append是添加的意思，就是把vc_current 添加（append） 到output_vcs这个空盒子里                                        
                                        # 就是将当前的爬升速度（climb_velocity）添加到 output_vcs 这个列表中
                                        # eg,在挑选每种水果（vc_current）之后，你将它们逐一放入购物袋（output_vcs 列表）
                                        # output_vcs.append(vc_current)： vc_current 就好比是一个水果（比如苹果），你首先把这个水果放入购物袋（output_vcs 列表）。
                                        
                                                                         
        if vc_current / vh >= 0:        # 然后，在购物袋中，你检查这个水果（vc_current）的属性，比如它是否新鲜（if vc_current / vh >= 0:）。            
            """Use vi/vh equation for positive vc"""
            """Power = T*(vc+vi)"""
            vi = -vc_current/2 + np.sqrt((vc_current/2)**2  +vh**2)
            P_act = Thrust*(vc_current + vi)/0.75
            vivh = vi/vh
            vc_vh = vc_current / vh
            
        elif vc_current / vh == 0:
            """Power = T*vh"""
            P_act = Thrust* vh/0.75
            
            
        elif -2 < vc_current / vh < 0:
            """Use k1 k2 k3 k4 equation"""
            vc_vh = vc_current / vh
            vi = vh * (k+ k1*vc_vh + k2 * vc_vh**2 + k3 * vc_vh**3 + k4 * vc_vh**4)
            vivh = vi/vh
            
            """Power = T*(vc+vi)"""
            P_act = Thrust*(vc_current + vi)/0.75

        elif vc_current / vh <= -2:
            """Use vi/vh equation for negative vc"""
            vi = -vc_current/2 - np.sqrt(vc_current**2/4 - vh**2)
            vivh = vi/vh
            vc_vh = vc_current / vh
            """Power = T*(vc+vi)"""
            P_act = Thrust*(vc_current + vi)/0.75
            
            
        else:
            raise RuntimeError("You messed up dude")
            """You did something wrong if this ever executes"""
            
        # Add final power to output_power
        output_power.append(P_act)  
        output_vis.append(vi)  
        output_vivh.append(vivh) 
        output_vcvh.append(vc_vh)
        output_vh.append(vh)
        # 把上面的记录簿output_power = []，记录上P_act （每次的测量值）       # eg，把上面算的P_act，放入output_power的购物袋
        
print("Powers:")
print(output_power)

print("\n Climb Velocities:")
print(output_vcs)
print('\n')

print("vi/vh:")
print(output_vivh)
print('\n')

print("vc/vh:")
print(output_vcvh)
print('\n')

print("vh:")
print(output_vh)
print('\n')



########################################################
conditions = ["Sea level", "5k feet", "10k feet"]

for i in range(0, len(conditions)):
    condition = conditions[i]
    start_index = i * 5
    vc_values = output_vcs[start_index:start_index + 5]
    power_values = output_power[start_index:start_index + 5]
    
    print(f"At {condition},")
    print(f"vc = {vc_values} ft/s,")
    print(f"Actual Power = {[round(power, 5) for power in power_values]} lbf*ft/s \n")
    
    
# plot:    
    
plt.plot(output_vcs[0:5], output_power[0:5])
plt.plot(output_vcs[5:10], output_power[5:10])
plt.plot(output_vcs[10:], output_power[10:])
plt.title("Power Required vs Climb Velocity")
plt.legend(["Sea level", "5k", "10k"])
plt.grid(True)
plt.xlabel("Vc (ft/s）")
plt.ylabel("Actual Power (lbf*ft/s)")
plt.show()


'''
这个代码作业强调了编程中常见的核心概念，包括变量初始化、循环、条件、列表操作和数学计算
初始化变量：了解如何初始化常量和变量是一个关键概念。在代码中，你初始化了给定的常量和临时变量，这些变量在后续计算中被使用。
循环结构：学习如何使用循环结构，包括嵌套循环。这个代码中使用了外部循环和内部循环，分别用于处理不同的条件和迭代。
列表操作：熟悉如何使用列表（output_rho、output_vcs、output_power）来存储和组织数据。你将不同的结果存储在列表中，以备将来使用


当整个代码被类比成一次超市购物经历，会更容易理解：

1.超市信息：首先，你来到了一个超市（代码的开头），上面写着“Problem 1”。这个超市卖的是直升机，你需要对它进行一些计算。

2.购物清单：你有一个购物清单（给定的常量和变量）：
    你知道直升机的重量（W）。
    你知道直升机可以以不同的爬升速度（Vc）飞行。
    你知道主旋翼的半径（R）和效率（FM）。
    你知道不同高度的空气密度（rho）。

3. 购物袋：你准备了三个购物袋（output_rho、output_vcs、output_power）来装购物。这些购物袋将用于存储不同的结果。

4. 购物开始：接下来，你开始购物：    
    你首先在不同的海拔条件下（密度变化）挑选水果（执行外部循环）。你挑选了海拔高度为海平面、5000英尺和10000英尺的水果。
    然后，你站在水果架前，挑选了不同的水果（执行内部循环）。这些水果代表了不同的爬升速度（Vc）。
    你把每个挑选的水果（vc_current）逐一放入购物袋（output_vcs 列表）。

5. 水果检查：然后，你开始检查每个水果的属性（执行条件语句）：
    如果水果是新鲜的（Vc/Vh >= 0），你使用一种特定的公式来计算功率，然后把这个水果（功率值）放入购物袋（output_power 列表）。
    如果水果已经检查过了（Vc/Vh == 0），你使用另一种公式来计算功率，然后也放入购物袋。
    如果水果有点特殊（-2 < Vc/Vh < 0），你使用一组复杂的公式，计算功率并放入购物袋。
    如果水果有问题（Vc/Vh <= -2），你使用另一种公式，计算功率，并也放入购物袋。
    如果你弄错了（其他情况），你会感到不安，但你抛出一个错误并继续购物。

6. 购物结果：最后，你完成了购物，看看购物袋中有什么：  
    你检查了购物袋中的功率（output_power 列表），看到了不同爬升速度下的功率需求。
    你还检查了购物袋中的水果清单（output_vcs 列表），这些是你挑选的不同爬升速度。
    最后，你展示了你的购物结果，就像在结账台上展示你买的水果一样（使用 plt 绘图）。
    
这个超市购物经历就是整个代码的执行过程，其中购物清单、购物袋、挑选水果和计算功率都是代码中的重要元素。希望这个比喻有助于更好地理解整个代码的逻辑。
'''

# Problem 2&3:

def fixed_point_iteration(mu_values, CT, alpha, max_iter=100, tol=1e-5):
    """
    Perform fixed-point iteration to calculate inflow ratios.

    Parameters:
        mu_values: Advanced ratio..
        CT: Thrust coefficient.
        alpha: Angle of attack in degrees.
        max_iter (int): Maximum number of iterations. Default is 100.
        tol (float): Tolerance level for convergence. Default is 1e-5.

    Returns:
        list: List of calculated inflow ratios.
    """
    inflow_ratios_fix_point = []

    for i in range(len(mu_values)):
        mu = mu_values[i]
        lam = np.sqrt(CT / 2)  # initialize lambda value

        for iteration in range(max_iter):
            new_lam = mu * np.tan(np.deg2rad(alpha)) + CT / (2 * np.sqrt(mu**2 + lam**2))

            if abs(new_lam - lam) < tol:    # 收敛足够                
                break

            lam = new_lam
            
        inflow_ratios_fix_point.append(lam)
        
    return inflow_ratios_fix_point


def newton_raphson_iteration(mu_values, CT, alpha, max_iter=100, tol=1e-5):

    inflow_ratios_newton_raphson = []

    for i in range(len(mu_values)):
        lam = np.sqrt(CT / 2)  # 初始化lambda值,lam_i
        mu = mu_values[i]

        for iteration in range(max_iter):
            f_lam = lam - (mu * np.tan(np.deg2rad(alpha))) - (CT / (2 * np.sqrt(mu**2 + lam**2)))
            g_lam = 1 + (CT*lam/2) * ((mu**2 + lam**2) ** (-3/2))
            new_lam = lam - f_lam / g_lam

            # 如果新旧值之间的差异小于容差，则认为达到了收敛
            if abs(new_lam - lam) < tol:
                break

            lam = new_lam

        inflow_ratios_newton_raphson.append(lam)
        
    return inflow_ratios_newton_raphson


# Usage:
CT = 0.01
mu_values = np.arange(0, 0.6, 0.05)
lamh = np.sqrt(CT / 2)
alphas = np.arange(-4, 9, 2)  # 攻角从-4度到8度，步长为2度

#################################################################output_vh
# Problem2: Plot for Fixed-point iteration and  Newton Raphson iteration when alpha =0
exact_result = np.sqrt(np.sqrt(1/4 * (mu_values / lamh) ** 4 + 1) - 1/2 * (mu_values / lamh) ** 2)
inflow_ratios_fp = fixed_point_iteration(mu_values, CT, alpha=0)
inflow_ratios_nr = newton_raphson_iteration(mu_values, CT, alpha=0)

plt.figure()
plt.plot(mu_values / lamh, exact_result, label='Exact Result', linestyle='-')
plt.scatter(mu_values / lamh, np.array(inflow_ratios_fp)/lamh, marker='o', label='Fixed-Point Method')
plt.scatter(mu_values / lamh, np.array(inflow_ratios_nr)/lamh, marker='*', label='Newton—Raphson Method')

plt.xlabel('mu/λ')
plt.ylabel('λi/λh')
plt.title('Inflow Ratios vs. Forward Speed Ratio for alpha = 0')
plt.legend()
plt.grid(True)
plt.show()

##################################################################
#Problem3:  Plot for Fixed-point iteration at different alpha
plt.figure()
plt.title('Fixed-point iteration: Inflow Ratios vs Mu')
plt.xlabel('Mu')
plt.ylabel('Inflow Ratios')
plt.grid(True)

plt.plot(mu_values / lamh, exact_result, label='Exact Result', linestyle='-')

for alpha in alphas:
    inflow_ratios_fx = fixed_point_iteration(mu_values, CT, alpha) 
    plt.plot(mu_values/lamh, inflow_ratios_fx/lamh, marker='*', linestyle='--', label=f'Alpha = {alpha}')
    
plt.legend()    # 添加图例
plt.show()

##################################################################
#Problem 3: Plot for Newton Raphson iteration at different alpha
plt.figure()
plt.title('Newton Raphson iteration: Inflow Ratios vs Mu')
plt.xlabel('Mu')
plt.ylabel('Inflow Ratios')
plt.grid(True)

plt.plot(mu_values / lamh, exact_result, label='Exact Result', linestyle='-')

for alpha in alphas:
    inflow_ratios_nr = newton_raphson_iteration(mu_values, CT, alpha)  # debug：参数顺序别错了
    plt.plot(mu_values/lamh, inflow_ratios_nr/lamh, marker='.', linestyle='--', label=f'Alpha = {alpha}')
    
plt.legend()
plt.show()

###################################################################
# Problem 4: Induced velicity curves 感应速度曲线

# print("vi/vh:")
# print(output_vivh)
# print('\n')

# print("vc/vh:")
# print(output_vcvh)
# print('\n')

# print("vh:")
# print(output_vh)
# print('\n')


# # 使用matplotlib绘制图像
# xx = (-4*vh, 4*vh, 400)
# # plt.plot(output_vcs/vh , output_vis/vh, 'r-', label='Induced Velocity Ratio')
# plt.plot(xx , output_vivh, 'b*', label='Induced Velocity Ratio')
# # plt.plot(np.array(output_vcs)/vh , np.array(output_vis)/vh, 'r-', label='Induced Velocity Ratio (Red Line)')
# plt.grid(True)
# plt.xlabel('Climb velocity ratio, Vc/Vh')
# plt.ylabel('Induced velocity ratio, Vi/Vh')
# plt.title('Induced velocity Curves')
# plt.legend()
# plt.show()



##
def calculate_vi(vc_current):
    vc_vh_ratio = vc_current / vh
    
    if vc_vh_ratio >= 0:
        return -vc_current/2 + np.sqrt((vc_current/2)**2 + vh**2)
    
    elif -2 < vc_vh_ratio < 0:
        return vh * (k + k1*vc_vh_ratio + k2*vc_vh_ratio**2 + k3*vc_vh_ratio**3 + k4*vc_vh_ratio**4)
    
    elif vc_vh_ratio <= -2:
        return -vc_current/2 - np.sqrt(vc_current**2/4 - vh**2)
    
    else:
        raise RuntimeError("Invalid vc/vh ratio")

# 创建x轴数据
vh=vh/35
vc = np.linspace(-4*vh, 4*vh, 1500)
# print('vh')
# print(vh)
# 使用上述函数来计算y轴数据 (感应速度)
vi = [calculate_vi(vc_val) for vc_val in vc]
# print(len(vc))
# print(len(vi))
# print("vi:")
# print(output_vis)

# # 使用matplotlib绘制图像
# plt.figure(figsize=(10,6))
# plt.plot(vc, vi, 'r-', label='complete induced velocity curve')

# plt.grid(True)
# plt.xlabel('Climb velocity ratio, Vc/Vh')
# plt.ylabel('Induced velocity ratio, Vi/Vh')
# plt.title('Induced velocity Curves')
# plt.legend()
# plt.show()

plt.figure(figsize=(10,6))

# x < -2
mask1 = vc <= -2
plt.plot(vc[mask1], np.array(vi)[mask1], 'r--', label='x <= -2')
# -2 <= x < 0
mask2 = (-2 < vc) & (vc <= 0)
plt.plot(vc[mask2], np.array(vi)[mask2], '-', label='-2 < x <= 0')

# x >= 0
mask3 = vc >= 0
plt.plot(vc[mask3], np.array(vi)[mask3], 'b--', label='x >= 0')

plt.grid(True)
plt.xlabel('Climb velocity ratio, Vc/Vh')
plt.ylabel('Induced velocity ratio, Vi/Vh')
plt.title('Induced velocity Curves')
plt.text(-0.5, 2, 'Descent ←', horizontalalignment='right')
plt.text(0.5, 2, '→ Climb', horizontalalignment='left')
plt.legend()
plt.show()

##

# velocity
x = np.linspace(-4, 3, 50)
v_normal = -0.5 * x + 0.5 * np.sqrt(x**2 + 4)
v_vortex = -1.125 * x + 0.974 - 1.372 * x**2 - 1.718 * x**3 - 0.655 * x**4
v_windmill = -0.5 * x - 0.5 * np.sqrt(x**2 - 4)

v_normal[(x < 0) | (x > 4)] = np.nan
v_vortex[(x < -2) | (x > 0)] = np.nan
v_windmill[(x < -4) | (x > -2)] = np.nan

plt.figure(figsize=(10,6))
plt.plot(x, v_normal, '--b', label='Normal working state')
plt.plot(x, v_vortex, '-r', label='Exact soln: Vortex ring state')
plt.plot(x, v_windmill, '--og', label='Windmill brake state')
plt.title('Power required as a function of climb and descent velocity')
plt.xlabel('Climb velocity ratio, Vc/Vh', color='#1C2833')
plt.ylabel('Power ratio, P/P_H', color='#1C2833')
plt.axhline(0, color='black',linewidth=0.5)
plt.axvline(0, color='black',linewidth=0.5)
plt.grid(color = 'gray', linestyle = '--', linewidth = 0.5)
plt.legend(loc='upper right')
plt.xlim([-3, 3])
plt.text(-0.5, 2, 'Descent ←', horizontalalignment='right')
plt.text(0.5, 2, '→ Climb', horizontalalignment='left')
# plt.axis('equal')
# plt.ylim([-4, 4])

plt.show()


# Power
x = np.linspace(-4, 3, 50)

power_normal = 0.5 * x + 0.5 * np.sqrt(x**2 + 4)
power_vortex = -0.125 * x + 0.974 - 1.372 * x**2 - 1.718 * x**3 - 0.655 * x**4
power_windmill = 0.5 * x - 0.5 * np.sqrt(x**2 - 4)

power_normal[(x < 0) | (x > 3)] = np.nan
power_vortex[(x < -2) | (x > 0)] = np.nan
power_windmill[(x < -4) | (x > -2)] = np.nan


plt.figure(figsize=(8,6))
plt.plot(x, power_normal, 'ob', label='Normal working state')
plt.plot(x, power_vortex, '-r', label='Exact soln: Vortex ring state')
plt.plot(x, power_windmill, 'og', label='Windmill brake state')
plt.title('Power required as a function of climb and descent velocity')
plt.xlabel('Climb velocity ratio, Vc/Vh', color='#1C2833')
plt.ylabel('Power ratio, P/P_H', color='#1C2833')
plt.axhline(0, color='black',linewidth=0.5)
plt.axvline(0, color='black',linewidth=0.5)
plt.grid(color = 'gray', linestyle = '--', linewidth = 0.5)
plt.legend(loc='upper left')
plt.xlim([-4, 4])
plt.ylim([-4, 4])
plt.text(-0.5, 2, 'Descent ←', horizontalalignment='right')
plt.text(0.5, 2, '→ Climb', horizontalalignment='left')
plt.show()





    
    
