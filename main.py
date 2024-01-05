import numpy as np
from matplotlib import pyplot as plt

## 输入气动参数
n = 10000
T0 = 1073
P0 = 1.382*101.325e3
m = 20
c0 = 238.855
gamma = 1.4
R = 287.08

## 输入几何参数
s = 184.1994e-3
l = 232.3242e-3
thickness_max = 33.57e-3
t_cl = 0.2e-3
h = (1773.94-1635.61)*1e-3
R_tip0 = 1635.61e-3
R_hub0 = 1074.5e-3
R_tip1 = 1773.94e-3
R_hub1 = 1072.98e-3
R_tip2 = 1852.46e-3
R_hub2 = 1072.98e-3
thickness_trailing = 2.247e-3
o = 93.9924e-3
D_m0 = R_tip0+R_hub0
D_m1 = R_tip1+R_hub1
D_m2 = R_tip2+R_hub2
A0 = np.pi * (R_tip0**2-R_hub0**2)
A1 = np.pi * (R_tip1**2-R_hub1**2)
A2 = np.pi * (R_tip2**2-R_hub2**2)
beta0 = np.deg2rad(30)
beta1 = np.deg2rad(70)
beta2 = np.deg2rad(30)
alpha0 = np.deg2rad(60)
alpha1 = np.deg2rad(30)
alpha2 = np.deg2rad(60)


## 定义函数
def X_function(lamda, gamma):
    return lamda / (np.sqrt(1 - (gamma-1)/(gamma+1)*lamda**2))

def B4_function(alpha2, beta2):
    return np.sin(alpha2) / np.sin(beta2)

def B3_function(A2, A1, beta2, beta1):
    return (A2*np.sin(beta2)) / (A1*np.sin(beta1))

def Y_funcion(lamda, theta, gamma):
    return lamda * (1 - (gamma-1)/(gamma+1)*(lamda/theta)**2)**(gamma/(gamma-1)) / (1 - (gamma-1)/(gamma+1)*lamda**2)

def lamda_function(C, Tc, gamma, R):
    return C / np.sqrt((2*gamma) / (gamma+1) * R * Tc)

## 计算静叶出口无量纲绝对速度系数
lamda_c0 = lamda_function(c0, T0, gamma, R)
phi = 0.9
e = 0.001
lamda_min = 0; lamda_max = 1
lamda_c1 = (lamda_max+lamda_min)/2
f1_function = lambda lamda: (Y_funcion(lamda_c0, 1, gamma) - B3_function(A1,A0,alpha1, alpha0) * Y_funcion(lamda, phi, gamma))


while abs(lamda_min - lamda_max)>e:
    if f1_function(lamda=lamda_min)*f1_function(lamda=lamda_c1)<0:
        lamda_min = lamda_min; lamda_max = lamda_c1
    else:
        lamda_min = lamda_c1; lamda_max = lamda_max
    
    lamda_c1 = (lamda_max+lamda_min)/2
print('静叶出口无量纲绝对速度系数=',lamda_c1)

## 计算静叶出口无量纲相对速度系数
lamda_min = 0; lamda_max = 1
lamda_w1 = (lamda_max+lamda_min)/2
f2_function = lambda lamda: (X_function(lamda_c1, gamma) - B4_function(beta1, alpha1)*X_function(lamda,gamma))
if f2_function(lamda_min)==0:
    lamda_w1 = lamda_min
elif f2_function(lamda_max)==0:
    lamda_w1 = lamda_max
elif f2_function(lamda_min)*f2_function(lamda_max)>0:
    print('重新输入边界点')
else:
    while abs(lamda_min - lamda_max)>e:
        if f2_function(lamda=lamda_min)*f2_function(lamda=lamda_w1)<0:
            lamda_min = lamda_min; lamda_max = lamda_w1
        else:
            lamda_min = lamda_w1; lamda_max = lamda_max
        
        lamda_w1 = (lamda_max+lamda_min)/2
print('静叶出口无量纲相对速度系数=', lamda_w1)

## 计算动叶出口无量纲相对速度系数
psi = 0.9
lamda_min = 0; lamda_max = 1
lamda_w2 = (lamda_max+lamda_min)/2
f3_function = lambda lamda: (Y_funcion(lamda_w1, 1, gamma) - B3_function(A2,A1,beta2, beta1) * Y_funcion(lamda, psi, gamma))
while abs(lamda_min - lamda_max)>e:
    if f3_function(lamda=lamda_min)*f3_function(lamda=lamda_w2)<0:
        lamda_min = lamda_min; lamda_max = lamda_w2
    else:
        lamda_min = lamda_w2; lamda_max = lamda_max
    
    lamda_w2 = (lamda_max+lamda_min)/2
print('动叶出口无量纲相对速度系数=', lamda_w2)

## 计算动叶出口无量纲绝对速度系数
lamda_min = 0; lamda_max = 1
lamda_c2 = (lamda_max+lamda_min)/2
f4_function = lambda lamda: (X_function(lamda_w2, gamma) - B4_function(alpha2, beta2)*X_function(lamda,gamma))
if f4_function(lamda_min)==0:
    lamda_c2 = lamda_min
elif f4_function(lamda_max)==0:
    lamda_c2 = lamda_max
elif f4_function(lamda_min)*f4_function(lamda_max)>0:
    print('重新输入边界点')
else:
    while abs(lamda_min - lamda_max)>e:
        if f4_function(lamda=lamda_min)*f4_function(lamda=lamda_c2)<0:
            lamda_min = lamda_min; lamda_max = lamda_c2
        else:
            lamda_min = lamda_c2; lamda_max = lamda_max
        
        lamda_c2 = (lamda_max+lamda_min)/2
print('动叶出口无量纲绝对速度系数=', lamda_c2)
lamda = np.linspace(0,1,50)
# X = X_function(lamda, gamma)
# plt.plot(lamda, X)
# plt.xlim(0,1)
# plt.ylim(0,1.2)
# plt.show()
# Y = Y_funcion(lamda, theta=0.8, gamma=gamma)
# plt.plot(lamda, Y)
# plt.xlim(0,1)
# plt.ylim(0,0.7)
# plt.show()

f1 = f4_function(lamda)
plt.plot(lamda, f1)
plt.show()
## step2


