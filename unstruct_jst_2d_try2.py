# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 10:32:26 2020

@author: Shalabh saxena
unstruct, 2d jst
"""
import numpy as np
import matplotlib.pyplot as plt

# range of phi_x and phi_y
Aphi_x = np.linspace(0,np.pi,num=50)
Aphi_y = np.linspace(0,np.pi,num=50)

# bigger range of sigma and ep2
nx  = 50
ny  = 40
sg  = np.linspace(0,1.2,num=nx)
#ep2 = np.linspace(0.35,0.354,num=ny)
ep2 = np.linspace(0.0,1,num=ny)

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

# there should be a nested loop to test for all possible combination
# value of phi_x and phi_y. This is achieved by creating this mesh grid
phi_x,phi_y = np.meshgrid(Aphi_x,Aphi_y)

# begin the loop; sgx=sgy=sg is an assumption here
k4 = 1.0/256.0
sqrt3    = np.sqrt(3)
sqrt3By2 = sqrt3/2.0


for i in range(nx):
    for j in range(ny):
        ep2j = ep2[j]
        ep4j = max(0,k4-ep2j)
        
        # Qu
        Qu_re  = -sqrt3*np.sin(phi_y/2)*np.sin(sqrt3By2*phi_x) +              \
        np.cos(phi_y/2)*np.cos(sqrt3By2*phi_x) - np.cos(phi_y)
        Qu_im  = np.sin(phi_y/2)*np.cos(sqrt3By2*phi_x) + np.sin(phi_y) +     \
        sqrt3*np.sin(sqrt3By2*phi_x)*np.cos(phi_y/2)
        Qu     = (2.0/3.0)*(Qu_re + 1j*Qu_im)
        
        # D2u
        D2u_re = -np.sin(phi_y/2)*np.sin(sqrt3By2*phi_x) +                    \
        sqrt3*np.cos(phi_y/2)*np.cos(sqrt3By2*phi_x) + np.cos(phi_y) -        \
        sqrt3 - 1
        D2u_im = sqrt3*np.sin(phi_y/2)*np.cos(sqrt3By2*phi_x) - np.sin(phi_y) \
        + np.sin(sqrt3By2*phi_x)*np.cos(phi_y/2)
        D2u    = (4.0/3.0)*ep2j*(D2u_re + 1j*D2u_im)

        # D4u
        D4u_re = 4*np.sin(phi_y/2)*np.sin(sqrt3By2*phi_x) -                   \
        np.sin(phi_y)*np.sin(sqrt3*phi_x) + 2*np.cos(phi_y)**2 -              \
        4*sqrt3*np.cos(phi_y/2)*np.cos(sqrt3By2*phi_x) +                      \
        sqrt3*np.cos(phi_y)*np.cos(sqrt3*phi_x) - 4*np.cos(phi_y) + 2 +3*sqrt3
        D4u_im = -(6*sqrt3+4)*np.sin(phi_y/2)*np.cos(sqrt3By2*phi_x) -         \
        2*np.sin(phi_y)*np.cos(phi_y) + sqrt3*np.sin(phi_y)*np.cos(sqrt3*phi_x)\
        + (2*sqrt3+4)*np.sin(phi_y) - 2*np.sin(sqrt3By2*phi_x)*np.cos(phi_y/2) \
        + np.sin(sqrt3*phi_x)*np.cos(phi_y)
        D4u    = (4.0/3.0)*ep4j*(D4u_re + 1j*D4u_im)
        
        # scaling of dissipation terms
        #D2u = D2u*2.857
        #D4u = D4u*2.857
        
        # alpha*dt
        alp = sg[i]*(Qu - D2u + D4u)
        
        # quantitities
        G    = 1 - alp + alp**2/2 - alp**3/4
        # RK3 is not giving any stable region, trying rk5
        #G    = 1 - alp + alp**2/2 - (3.0/16.0)*alp**3 + (3.0/96.0)*alp**4 -   \
        #(3.0/384.0)*alp**5
        # Euler step
        #G    = 1 - alp
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
levels = [0.0,1.0]
#levels = np.linspace(0,1.0,num=20)
cp = plt.contourf(xx,yy,Edmax,cmap='hsv',levels=levels)
#cp = plt.contourf(xx,yy,Edmax,cmap='hsv')
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
    

# plot |G| at ep2=0.2 sq=0.5
ep2j = 0.2
ep4j = max(0,k4-ep2j)

# Qu
Qu_re  = -sqrt3*np.sin(phi_y/2)*np.sin(sqrt3By2*phi_x) +              \
np.cos(phi_y/2)*np.cos(sqrt3By2*phi_x) - np.cos(phi_y)
Qu_im  = np.sin(phi_y/2)*np.cos(sqrt3By2*phi_x) + np.sin(phi_y) +     \
sqrt3*np.sin(sqrt3By2*phi_x)*np.cos(phi_y/2)
Qu     = (2.0/3.0)*(Qu_re + 1j*Qu_im)

# D2u
D2u_re = -np.sin(phi_y/2)*np.sin(sqrt3By2*phi_x) +                    \
sqrt3*np.cos(phi_y/2)*np.cos(sqrt3By2*phi_x) + np.cos(phi_y) -        \
sqrt3 - 1
D2u_im = sqrt3*np.sin(phi_y/2)*np.cos(sqrt3By2*phi_x) - np.sin(phi_y) \
+ np.sin(sqrt3By2*phi_x)*np.cos(phi_y/2)
D2u    = (4.0/3.0)*ep2j*(D2u_re + 1j*D2u_im)

# D4u
D4u_re = 4*sqrt3*np.sin(phi_y/2)*np.sin(sqrt3By2*phi_x) -             \
sqrt3*np.sin(phi_y)*np.sin(sqrt3*phi_x) - 2*np.cos(phi_y)**2 -        \
4*np.cos(phi_y/2)*np.cos(sqrt3By2*phi_x) +                            \
np.cos(phi_y)*np.cos(sqrt3*phi_x) + 4*np.cos(phi_y) + 1
D4u_im = -2*np.sin(phi_y/2)*np.cos(sqrt3By2*phi_x) + np.sin(2*phi_y)  \
+ np.sin(phi_y)*np.cos(sqrt3*phi_x) - 2*np.sin(phi_y) -               \
2*sqrt3*np.sin(sqrt3By2*phi_x)*np.cos(phi_y/2) +                      \
sqrt3*np.sin(sqrt3*phi_x)*np.cos(phi_y)
D4u    = (4.0/3.0)*ep4j*(D4u_re + 1j*D4u_im)

alp = 0.5*(Qu - D2u + D4u)
G    = 1 - alp + alp**2/2 - alp**3/4
Gmod = np.absolute(G)

plt.figure()
plt.xlabel('$\phi_x$')
plt.ylabel('$\phi_y$')
plt.grid(True)
levels=np.linspace(0.0,1.0,num=10)
cp = plt.contourf(phi_x,phi_y,Gmod,cmap='hsv',levels=levels)
#cp = plt.contourf(xx,yy,Edmax,cmap='hsv')
plt.colorbar(cp)
plt.title('Contours of |G| @ $\epsilon_2$ = 0.2 and $\sigma$ = 0.5')