# -*- coding: utf-8 -*-
"""
Created on Sun Jan 26 08:08:05 2020

@author: Shalabh saxena
structured,2d,full JST,contours
"""
import numpy as np
import matplotlib.pyplot as plt

# range of phi_x and phi_y
Aphi_x = np.linspace(0,np.pi,num=50)
Aphi_y = np.linspace(0,np.pi,num=50)

# instead of a nested loop over phi_x and phi_y value, create meshgrid
phi_x, phi_y = np.meshgrid(Aphi_x,Aphi_y)

# bigger range of sigma and ep2
nx  = 50
ny  = 40
sg  = np.linspace(0,1.2,num=nx)
ep2 = np.linspace(0,0.25,num=ny)
#ep2 = np.linspace(0,0.004,num=ny)  # only in 4th order dissipation region

# define max diffusion error
Edmax = np.zeros([nx,ny])

# Gmod avg
Agmod = np.zeros([nx,ny])

# define min diffusion error
Edmin = np.zeros([nx,ny])

# dispersion error
Edisp = np.zeros([nx,ny])

# create 2d mesh
xx,yy = np.meshgrid(sg,ep2,indexing='ij')

# begin the loop
k4 = 1.0/256.0
for i in range(nx):
    for j in range(ny):
        ep2j = ep2[j]
        ep4j = max(0,k4-ep2j)
        lx   = np.cos(2*phi_x) - 4*np.cos(phi_x) + 3
        ly   = np.cos(2*phi_y) - 4*np.cos(phi_y) + 3
        g1x  = ep2j*(1-np.cos(phi_x)) + ep4j*lx
        g1y  = ep2j*(1-np.cos(phi_y)) + ep4j*ly
        g2x  = 1j*np.sin(phi_x) + 2*g1x
        g2y  = 1j*np.sin(phi_y) + 2*g1y
        alp  = sg[i]*(g2x+g2y)
        G    = 1 - alp + alp**2/2 - alp**3/4
        GC   = np.conjugate(G)
        PHI  = np.angle(GC)/(sg[i]*(phi_x+phi_y)+0.000001)
        Gmod = np.absolute(G)
        Edmax[i,j] = np.max(Gmod)
        Edmin[i,j] = np.min(Gmod)
        Edisp[i,j] = np.max(abs(1-PHI))
        Agmod[i,j] = np.mean(Gmod)
    

#plot the results - max|G|
plt.figure()
plt.xlabel('Courant No.')
plt.ylabel('$\epsilon_2$')
plt.grid(True)
levels = [0.0,1]
cp = plt.contourf(xx,yy,Edmax,cmap='hsv',levels=levels)
plt.colorbar(cp)
plt.title('Contours of max |G|')

#plot the results - min|G|
plt.figure()
plt.xlabel('Courant No.')
plt.ylabel('$\epsilon_2$')
plt.grid(True)
levels = np.linspace(0.0,1.0,num=10) 
cp = plt.contourf(xx,yy,Edmin,cmap='hsv',levels=levels)
plt.colorbar(cp)
plt.title('Contours of min |G|')

#plot the results - max phase errors
plt.figure()
plt.xlabel('Courant No.')
plt.ylabel('$\epsilon_2$')
plt.grid(True)
cp = plt.contourf(xx,yy,Edisp,cmap='hsv')
plt.colorbar(cp)
plt.title('Contours of max phase errors')

# plot the average dissipation contours
plt.figure()
plt.xlabel('Courant No.')
plt.ylabel('$\epsilon_2$')
plt.grid(True)
levels = np.linspace(0.0,1.0,num=10)
cp = plt.contourf(xx,yy,Agmod,cmap='hsv',levels=levels)
plt.colorbar(cp)
plt.title('Contours of average |G|')