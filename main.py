import numpy as np
from matplotlib import pyplot as plt

## input dynamic parameters
n = 3600
T0 = 1633
P0 = 1.3e6
m = 26
c0 = 100
gamma = 1.4
R = 8.314

## input geometric parameters
s = 21
l = 20
thickness_max = 1.12
t_cl = 0.02
h = 18
R_tip = 273
R_hub = 197
thickness_trailing = 0.26
o = 8.11

## step1
def X_function(lamda, gamma):
    return lamda / (np.sqrt(1 - (gamma-1)/(gamma+1)*lamda**2))

def B4_function(alpha2, beta2):
    return np.sin(alpha2) / np.sin(beta2)

def B3_function(A1, A2, beta1, beta2):
    return (A2*np.sin(beta2)) / (A1*np.sin(beta1))

def Y_funcion(lamda, theta, gamma):
    return lamda * (1 - (gamma-1)/(gamma+1)*(lamda/theta)**2)**(gamma/(gamma-1)) / (1 - (gamma-1)/(gamma+1)*lamda**2)

def lamda_function(C, Tc, gamma, R):
    return C / np.sqrt((2*gamma) / (gamma+1) * (R*Tc))

lamda_c0 = lamda_function(c0, T0, gamma, R)
phi = 0
## step2


