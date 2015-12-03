# -*- coding: utf-8 -*-
"""
Created on Mon Oct 12 03:08:51 2015

@author: ilovealltigers
"""
from numpy import pi
import numpy as np
import scipy.constants as spc
import matplotlib.pylab as plt

globalvar = rho = np.float128(raw_input('Input Stellar Density in kilograms per cubic meter:'))

def main():
    RungeKutta4(f1, f2, 1e-15, 1000 , 0.0, 1.6e35, 1000)
    

def RungeKutta4(f, g, x0, x1, y0, z0, n):
    h = (x1-x0)
    x = x0
    y = y0
    z = z0
    data_y = []
    data_z = []
    #data=[]
    i = 0
    for i in range(0, n + 1):

        k1 = h * f(x,y)
        j1 = h * g(x,z)

        k2 = h * f(x + h/2.0, y + k1/2.0)
        j2 = h * g(x + h/2.0, z + j1/2.0)

        k3 = h * f(x + h/2.0, y + k2/2.0)
        j3 = h * g(x + h/2.0, z + j2/2.0)

        k4 = h * f(x + h, y + k3)
        j4 = h * g(x + h, z + j3)
        
        k = (1.0/6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4)
        j = (1.0/6.0) * (j1 + 2.0*j2 + 2.0*j3 + j4)   
        x = x + h
        y = y + k
        z = z + j
        #data.append([y,z])
        data_y.append(y)
        data_z.append(z)
    #print data
    plt.plot(data_y,data_z)
    plt.xlabel('Mass (kilograms)', fontsize=14)
    plt.ylabel('Pressure (pascals) ', fontsize=14)
    plt.show()    

def f1(r,m):
    return 4.0 * pi * (r**2) * rho
    
def f2(r,P):
    return ((-4.0/3.0) * spc.G * r * rho**2)  


if __name__ == '__main__':
   main()

