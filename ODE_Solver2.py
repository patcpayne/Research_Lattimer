# -*- coding: utf-8 -*-
"""
Created on Mon Oct 12 03:08:51 2015

@author: Patrick Payne

Purpose: Simple Hydrostatic Equilibrium: Star Model
"""

from numpy import pi
import decimal as d
import numpy as np
#import scipy.constants as spc
import matplotlib.pylab as plt

#globalvar = rho = np.float128(raw_input('Input Stellar Density (kilograms per cubic meter):'))
#rho = rho * (spc.G/(spc.c**2)) * 1.0e6

def main():
    RungeKutta4(f1, f2, 0.0, 1000, 0.0, 2.041e-22, 1.36, 10000)
    RK4(1.36, 10000)

"""-------------------------------------------------------------------------"""

def RungeKutta4(f, g, x0, x1, y0, z0, gamma, n):
    h = x1  # Step Size
    x = x0  # Initial Radius
    y = y0  # Initial Mass 
    z = z0  # Central Pressure 
    rho = z0 ** (1.0/gamma)
    data_x = []
    data_y = []
    data_z = []
    data=[]
    i = 0
    for i in range(0, n + 1):
        if x < 1.0e-25:
            x = x + 1.0e-10
            y = y
            z = z
        else:    
            k1 = h * f(x,y, rho)
            k2 = h * f(x + h/2.0, y + k1/2.0, rho)
            k3 = h * f(x + h/2.0, y + k2/2.0, rho)            
            k4 = h * f(x + h, y + k3, rho)            
            k = (1.0/6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4)
            
            j1 = h * g(x, z, rho)
            j2 = h * g(x + h/2.0, z + j1/2.0, rho)
            j3 = h * g(x + h/2.0, z + j2/2.0, rho)
            j4 = h * g(x + h, z + j3, rho)
            j = (1.0/6.0) * (j1 + 2.0*j2 + 2.0*j3 + j4)   
            print rho, x, y, z+j
            if (z + j) > 0.0:                    # Determination if step size needs
                x = x + h                        # to change for next iteration
                y = y + k                        
                z = z + j
                rho = z ** (1.0/gamma)
            elif (z + j) < 0.0:
                h = (-1.0) * (z / j) * h
                x = x
                y = y
                z = z
            if z < 1.0e-25:
                x = x + h
                y = y + k
                z = z + j
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
    plt.show() 
    x = x * 1.0
    y = y * 1.0
    print "Radius (km)= %g" %x
    print "Mass (km)= %.4E" %d.Decimal(str(y)) #Converts to Scientific Notation
    

def f1(r, m, rho): 
    dM_dr = 4.0 * pi * (r**2) * rho# Function for dM/dr
    return   dM_dr
  
def f2(r, P, rho):            # Function for dP/dr, substituted m with Integral(dm/dr) 
    m =  (4.0/3.0) * pi * r**3 * rho      
    dP_dr = (-1.0) * ( m * rho) / r**2
    return  dP_dr   
    
'''-------------------------------------------------------------------------'''

def Initials(gamma):
    r = 0.0
    dr = 1000.0
    M0 = 0.0
    P0 = 2.041e-22
    Rho0 = P0**(1.0/gamma)
    y0 = np.array([M0, P0])
    return (r, dr, y0, Rho0)
    
def Eqns(r,y0,Rho):
    M = y0[0]
    P = y0[1]
    rhs1 = 4.0 * pi * (r**2) * Rho 
    rhs2 = (-1.0) * ( (4.0/3.0) * pi * r**3 * Rho * Rho) / r**2 
    rhs = np.array([rhs1, rhs2])    
    return  M, P, rhs 
    
def RK4(gamma, n):
    r, dr, y0 , rho = Initials(gamma)
    data_r = []
    data_M = []
    data_P = []
    y = y0

    for i in range(0, n + 1):
        
        if r < 1.0e-25:
            r = r + 1.0e-20
            y = y
        
        M , P , rhs = Eqns(r,y0,rho)
        k1 = dr * rhs
        
        M , P , rhs =Eqns(r + dr/2.0, y0 + k1/2.0, rho)
        k2 = dr * rhs
        
        M , P , rhs =Eqns(r + dr/2.0, y0 + k2/2.0, rho)
        k3 = dr * rhs
        
        M , P , rhs =Eqns(r + dr, y + k3, rho)
        k4 = dr * rhs
        
        k = (1.0/6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
        j = k[1]        
        
        if (P + j) > 0.0:
            y = y + k
            r = r + dr
            rho = y[1]**(1.0/gamma)
            data_r.append(r)
            data_M.append(M)
            data_P.append(P)
        elif (P + j) < 0.0:
            dr = (P / j) * dr * (-1.0)
            y = y
        if P < 1.0e-25:
            r = r + dr
            y = y + k 
            data_r.append(r)
            data_M.append(M)
            data_P.append(P)
            break   
    plt.plot(data_r,data_M)                  # Plots Radius v Mass 
    plt.xlabel('Radius ($km$)', fontsize=14)                     
    plt.ylabel('Mass ($km$) ', fontsize=14)
    plt.show()
    plt.plot(data_r,data_P)                  # Plots Radius v Pressure
    plt.xlabel('Radius ($km$)', fontsize=14)                     
    plt.ylabel('Pressure ($km^{-2}$) ', fontsize=14)    
    plt.show()
    plt.plot(data_M,data_P)                  # Plots data (Mass v Pressure)
    plt.xlabel('Mass ($km$)', fontsize=14)                     
    plt.ylabel('Pressure ($km^{-2}$) ', fontsize=14)
    plt.show() 
    print "Radius (km)= %g" %r
    print "Mass (km)= %.4E" %d.Decimal(str(M)) #Converts to Scientific Notation
            
"""-------------------------------------------------------------------------"""            

if __name__ == '__main__':
   main()# Research_Lattimer