# -*- coding: utf-8 -*-
"""
Created on Mon Oct 12 03:08:51 2015

@author: Patrick Payne

Purpose: Simple Hydrostatic Equilibrium: Star Model
"""

from numpy import pi
import decimal as d
import numpy as np
import scipy.constants as spc
import matplotlib.pylab as plt

globalvar = rho = np.float128(raw_input('Input Stellar Density (kilograms per cubic meter):'))
rho = rho * (spc.G/(spc.c**2)) * 1.0e6

def main():
    RungeKutta4(f1, f2, 1.0e-163, 1000.0, 0.0, 2.041e-22, 1000)

def RungeKutta4(f, g, x0, x1, y0, z0, n):
    h = x1  # Step Size
    x = x0  # Initial Radius
    y = y0  # Initial Mass 
    z = z0  # Central Pressure  
    data_x = []
    data_y = []
    data_z = []
    data=[]
    i = 0
    for i in range(0, n + 1):
        if x < 1.0e-25:
            x = x + h
            y = y
            z = z
        else:    
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
           #5print x, y, z
            if (z + ((1.0/6.0) * (j1 + 2.0*j2 + 2.0*j3 + j4))) > 0.0:                    # Determination if step size needs
                x = x + h                        # to change for next iteration
                y = y + (1.0/6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4)                        
                z = z + (1.0/6.0) * (j1 + 2.0*j2 + 2.0*j3 + j4) 
            elif (z + ((1.0/6.0) * (j1 + 2.0*j2 + 2.0*j3 + j4))) < 0.0:
                h = (-1.0) * (z / ((1.0/6.0) * (j1 + 2.0*j2 + 2.0*j3 + j4)) ) * h
                x = x
                y = y
                z = z
            if z < 1.0e-25:
                x = x + h
                y = y + (1.0/6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4)
                z = z + (1.0/6.0) * (j1 + 2.0*j2 + 2.0*j3 + j4) 
                break
        data.append([x,y,z])
        data_x.append(x)
        data_y.append(y)
        data_z.append(z)
      
# The units for Mass and Pressure are geometric units obtained by
 #      multiplying MASS by G/(c**2) and PRESSURE by G/(c**4) and then 
  #     converting from meters to kilometers dividing the by 1000 for MASS and
   #    multiplying by 1000000 for PRESSURE
    #   Gravitational constant and Speed of light drop out of equations as they
     #  equal 1 in the geometrized unit system 
       
    plt.plot(data_x,data_y)                  # Plots Radius v Mass 
    plt.xlabel('Radius ($km$)', fontsize=14)                     
    plt.ylabel('Mass ($km$) ', fontsize=14)
    plt.show()
    plt.plot(data_x,data_z)                  # Plots Radius v Pressure
    plt.xlabel('Radius ($km$)', fontsize=14)                     
    plt.ylabel('Pressure ($km^{-2}$) ', fontsize=14)    
    plt.show()
    plt.plot(data_y,data_z)                  # Plots data (Mass v Pressure)
    plt.xlabel('Mass ($km$)', fontsize=14)                     
    plt.ylabel('Pressure ($km^{-2}$) ', fontsize=14)
    plt.plot(0)
    plt.show() 
    x = x * 1.0
    y = y * 1.0
    print "Radius (km)= %g" %x
    print "Mass (km)= %.2E" %d.Decimal(str(y)) #Converts to Scientific Notation

 
def f1(r,m):                                 # Function for dM/dr
    return   4.0 * pi * (r**2) * rho
  
def f2(r,P):           # Function for dP/dr, substituted m with Integral(dm/dr) 
    m =  (4.0/3.0) * pi * r**3 * rho         
    return   (-1.0) * ( m * rho) / r**2  
    
if __name__ == '__main__':
   main()# Research_Lattimer
