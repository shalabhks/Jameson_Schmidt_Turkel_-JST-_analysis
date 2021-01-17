# -*- coding: utf-8 -*-
"""
Created on Tue Apr  7 12:09:05 2020

@author: Shalabh saxena
"""

import numpy as np
import matplotlib.pyplot as plt

nx = 50
ny = 50 # lines
ncx = nx - 1
ncy = ny - 1
U   = np.zeros((ncx+4,ncy+4)) # 2 layers of aux cells
A   = 1.0
B   = 1.0
xd  = 2.0  # domain size in x from -1 to 1
yd  = 2.0  # domain size in y from -1 to 1
dx  = xd / (nx-1)
dy  = yd / (ny-1)
k2  = 1/4
k4  = 1/256
currentTime = 0.0
delu  = np.zeros((ncx+4,ncy+4))

# x distribution and y distribution
x = np.zeros(nx-1)
y = np.zeros(ny-1)
for ii in range(2,nx+1):
    i    = ii - 2
    x[i] = -1 + (2*i+1)/(nx-1)
for jj in range(2,ny+1):
    j    = jj - 2
    y[j] = -1 + (2*j+1)/(ny-1)

# initialize solution
for ii in range(2,nx+1):
    for jj in range(2,ny+1):
        i = ii - 2
        j = jj - 2
        r = np.sqrt(x[i]**2+y[j]**2)
        t = (r/0.2)*(np.pi/2)
        U[ii,jj] = np.cos(t)
        if r > 0.2: U[ii,jj] = 0.0
    
def contour_plot():
    """
    Plots contour of solution
    """
    xx,yy = np.meshgrid(x,y)
    plt.figure()
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.grid(True)
    levels = np.linspace(0,1.0,num=10)
    #cp = plt.contourf(xx,yy,U[2:nx+1,2:ny+1],cmap='hsv',levels=levels)
    cp = plt.contour(xx,yy,U[2:nx+1,2:ny+1],cmap='hsv',levels=levels)
    plt.clabel(cp,inline=True,fontsize=8)
    #plt.colorbar(cp)
    plt.title('Contour plot')

def plot_line(write):
    """
    line plot
    """
    z=np.zeros(nx-1)
    s=np.zeros(nx-1)
    for ii in range(2,nx+1):
        i = ii - 2
        z[i]=U[ii,ii]
        s[i]=np.sqrt((x[i]+1)**2+(y[i]+1)**2)
    str1 = 't = ' + str(currentTime)
    plt.plot(s,z,label=str1)
    plt.xlabel('S')
    plt.ylabel('U')
    plt.title('U along y=x line')
    plt.grid(True)
    plt.legend()
    
    # if write the write data to text file for comparison to destruct and 
    # unstruct
    if write:
        s = s.reshape(nx-1,1)
        z = z.reshape(nx-1,1)
        ar = np.hstack((s,z))
        np.savetxt('destruct.dat',ar)
    

# plot initial solution
fig_l = plt.figure(1)
plot_line(write=False)
contour_plot()
#

# calculate convective flux
def calc_qu():
    q = np.zeros((nx-1,ny-1))
    for ii in range(2,nx+1):
        for jj in range(2,ny+1):
            i = ii - 2
            j = jj - 2
            ue = 0.5*(U[ii+1,jj] + U[ii,jj])
            uw = 0.5*(U[ii-1,jj] + U[ii,jj])
            un = 0.5*(U[ii,jj+1] + U[ii,jj])
            us = 0.5*(U[ii,jj-1] + U[ii,jj])
            q[i,j] = A*ue*dy -A*uw*dy + B*un*dx - B*us*dx
        
    return q

# calc_d2u
def calc_d2u():
    d2u = np.zeros((nx-1,ny-1))
    for ii in range(2,nx+1):
        for jj in range(2,ny+1):
            i  = ii - 2
            j  = jj - 2
            nui = np.abs(U[ii+1,jj]) + np.abs(U[ii-1,jj]) - 2*np.abs(U[ii,jj])
            dei = np.abs(U[ii+1,jj]) + np.abs(U[ii-1,jj]) + 2*np.abs(U[ii,jj])
            nuj = np.abs(U[ii,jj+1]) + np.abs(U[ii,jj-1]) - 2*np.abs(U[ii,jj])
            dej = np.abs(U[ii,jj+1]) + np.abs(U[ii,jj-1]) + 2*np.abs(U[ii,jj])
            nu  = nui + nuj
            de  = dei + dej + 1.0e-6
            mu  = np.abs(nu)/de
            due = U[ii+1,jj] - U[ii,jj]
            duw = U[ii-1,jj] - U[ii,jj]
            dun = U[ii,jj+1] - U[ii,jj]
            dus = U[ii,jj-1] - U[ii,jj]
            d2u[i,j] = mu*(A*dy*due + A*dy*duw) + mu*(B*dx*dus + B*dx*dun)
            d2u[i,j] = k2*d2u[i,j]
            delu[ii,jj] = due + duw + dun + dus
        
    return d2u

def calc_d4u():
    d4u = np.zeros((nx-1,ny-1))
    
    for ii in range(2,nx+1):
        for jj in range(2,ny+1):
            i  = ii - 2
            j  = jj - 2
            nu = np.abs(U[ii+1,jj]) + np.abs(U[ii-1,jj]) - 2*np.abs(U[ii,jj])
            de = np.abs(U[ii+1,jj]) + np.abs(U[ii-1,jj]) + 2*np.abs(U[ii,jj])
            nu = nu + np.abs(U[ii,jj+1]) + np.abs(U[ii,jj-1]) - 2*np.abs(U[ii,jj])
            de = de + np.abs(U[ii,jj+1]) + np.abs(U[ii,jj-1]) + 2*np.abs(U[ii,jj])
            mu = np.abs(nu)/(de+1.0e-6)
            ep4 = max(0,(k4-mu*k2))
            due  = delu[ii+1,jj] - delu[ii,jj] 
            duw  = delu[ii-i,jj] - delu[ii,jj]
            dun  = delu[ii,jj+1] - delu[ii,jj]
            dus  = delu[ii,jj-1] - delu[ii,jj]
            d4u[i,j] = ep4*(A*dy*(due + duw) + B*dx*(dun+dus)) 
       
    return d4u
            
# apply periodic BC
def apply_bc():
    U[0,: ] = U[-4,:]
    U[1,: ] = U[-3,:]
    U[-2,:] = U[2,: ]
    U[-1,:] = U[3,: ]
    U[:,0 ] = U[:,-4]
    U[:,1 ] = U[:,-3]
    U[:,-2] = U[:,2 ]
    U[:,-1] = U[:,3 ]
    
    delu[0,: ] = delu[-4,:]
    delu[1,: ] = delu[-3,:]
    delu[-2,:] = delu[2,: ]
    delu[-1,:] = delu[3,: ]
    delu[:,0 ] = delu[:,-4]
    delu[:,1 ] = delu[:,-3]
    delu[:,-2] = delu[:,2 ]
    delu[:,-1] = delu[:,3 ]
    

# residue plot
plt.figure(3)
plt.title("Residue Vs Iter")
plt.xlabel('Iterations')
plt.ylabel('L2 Norm of del U')

# march in time
istep = 0
nstep = 2000
dt    = 1.0e-3
area  = (2*2)/((nx-1)*(ny-1))
R0    = np.zeros((nx-1,ny-1))

while istep < nstep:
    istep = istep + 1
    
    # Runge Kutta steps
    Ucopy = U.copy()

    for ik in range(3):
        Qu    = calc_qu()
        d2u   = calc_d2u()
        d4u   = calc_d4u()
        R     = (-Qu + d2u - d4u)/(dx*dy)
        
        if ik==0:
            R0 = R.copy()
        if ik==0:
            U[2:nx+1,2:ny+1] = Ucopy[2:nx+1,2:ny+1] + dt*R
        elif ik==1 or ik==2:
            U[2:nx+1,2:ny+1] = Ucopy[2:nx+1,2:ny+1] + (dt/2)*(R + R0)
        else:
            print('Only RK3 implemented!')
        
        apply_bc()
        
    res = np.sum(R**2)
    res = np.sqrt(res/(nx*ny))
    
    print("Iter = ",istep," Residue = ",res)
    currentTime = currentTime + dt
    plt.figure(3)
    plt.plot(istep,res,'b.')
    plt.pause(0.05)
    

plt.figure(1)
plot_line(write=True)
#contour_plot()


    
        
    