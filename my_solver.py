# -*- coding: utf-8 -*-
"""
Created on Sun Mar 29 11:19:51 2020

@author: Shalabh saxena
"""
import numpy as np
#from scipy.interpolate import Rbf
from scipy import interpolate
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

class cell_ds:
    value   = 0.0
    center  = np.zeros(2)
    triIndx = 0
    aux     = False
    area    = 0.0
    
#class b2b_ds():
#    j0 = 0    # aux cell
#    j1 = 0    # interior cell
#    jn = 0    # connected cell to get value from

class solver:
    
    U   = []  # Cell to store Scalar that is transported
#    B2B = []  # B2B connectivity for periodic BC
    nn     = 0
    nn_int = 0
    nptsNew = 0
    A      = 0 # X- Convective velocity
    B      = 0 # Y- Convective Velocity
    cf     = np.array([]) # convective flux
    d2f    = np.array([]) # 2nd order diffusive flux
    k2     = 1/4         # k2
    k4     = 1/256
    Ctime  = 0.0          # current time of iteration
    cf     = np.array([]) # convective flux
    d2f    = np.array([]) # diffusive flux 2nd order
    d4f    = np.array([]) # 4th order diffusive flux
    sw     = np.array([]) # switch
    scale  = 0.0          # line plot is not 1.0 in initial data so scaling it
                          # up for visualization only
    flagWrite = False     # If true line plot writes out a text file for comparative plots with
                          # structured and destructured
    
    def __init__(self,mesh):
        self.nn_int   = mesh.nelem
        self.nn       = mesh.nelem + mesh.nBCedges
        self.nptsNew  = mesh.nelem
        self.A        = 1.0
        self.B        = 1.0
        self.Ctime    = 0.0
        #rmax          = np.sqrt(mesh.xmax**2 + mesh.ymax**2)
        
        # create interior cells
        i = 0
        for t in mesh.tri:
            u         = cell_ds()
            u.value   = 0.0
            u.triIndx = i
            u.center  = (1/3)*(mesh.points[t[0]]+mesh.points[t[1]]+mesh.points[t[2]])
            xc        = u.center[0]
            yc        = u.center[1]
            rc        = np.sqrt(xc**2 + yc**2)
            theta     = (rc/0.2)*(np.pi/2)
            u.value   = np.cos(theta)
            u.aux     = False
            if rc > 0.2: 
                u.value = 0.0
            v  = mesh.points[t[1]]-mesh.points[t[0]]
            s1 = np.sqrt(np.inner(v,v))
            v  = mesh.points[t[2]]-mesh.points[t[1]]
            s2 = np.sqrt(np.inner(v,v))
            v  = mesh.points[t[2]]-mesh.points[t[0]]
            s3 = np.sqrt(np.inner(v,v))
            s  = (s1+s2+s3)/2
            u.area = np.sqrt(s*(s-s1)*(s-s2)*(s-s3))
                                    
            self.U.append(u)
            i         = i + 1
            
        # create aux cells and b2b index
#        for e in mesh.edges:
#            if e.bndry == 0:
#                continue
#            u         = cell_ds()
#            u.center  = np.zeros(2)
#            u.triIndx = i
#            i         = i + 1
#            u.value   = 0.0
#            u.aux     = True
#            u.area    = 0.0
#            self.U.append(u)
            
#            b2b    = b2b_ds()
#            b2b.j1 = e.t1
#            b2b.jn = -e.t2
#            b2b.j0 = i-1
#            self.B2B.append(b2b)
            
    # plot line graph at 45 degrees
    def plot_line_graph(self,mesh,iFig):
        """
        Plots the distribution of U along line y=x
        """
        # extract information
        xb = np.zeros((1,2))
        zb = np.array([])
        for u in self.U:
            if u.aux:
                continue
            xb[-1] = u.center
            xb     = np.vstack((xb,np.zeros((1,2))))
            zb = np.append(zb,u.value)
        xb = xb[:-1,:]
        
        # build points to interpolate data on
        x = np.linspace(mesh.xmin,mesh.xmax,num=50)
        x = x.reshape(x.shape[0],1)
        x = np.hstack((x,x))
        s = np.sqrt((x[:,0]-mesh.xmin)**2+(x[:,1]-mesh.ymin)**2)
        
        # interpolate data 1D
        zb = zb.reshape(zb.shape[0],1)
        f = interpolate.LinearNDInterpolator(xb,zb)
        z = f(x)
        
        # get scaling factor; Interpolation on initial data also does not
        # capture the peak (1.0), so for comparison with struct etc we ned
        # to scale the data by this cons factor
        if self.scale == 0:
            self.scale = 1.0/np.max(z[1:-2])
            print ("Scale factor =",self.scale)
        
        # scale values
        z = z * self.scale
        
        # plot the figure
        plt.figure(iFig)
        plt.grid(True)
        str1 = "t = " + str(self.Ctime)
        plt.plot(s,z,label=str1)
        plt.xlabel('s')
        plt.ylabel('U')
        plt.title('U along y=x line')
        plt.legend()
        
        # write out file if asked
        if self.flagWrite:
            s = s.reshape(s.shape[0],1)
            z = z.reshape(z.shape[0],1)
            a = np.hstack((s,z))
            np.savetxt('unstruct.dat',a)
        
           
    # plot contour of U    
    def plot_contour(self,mesh):
        """
        plots the contour of distribution on U
        """
        # get data
        xb = np.array([])
        yb = np.array([])
        zb = np.array([])
        for u in self.U:
            if u.aux:
                continue
            xb = np.append(xb,u.center[0])
            yb = np.append(yb,u.center[1])
            zb = np.append(zb,u.value)
        #zb = np.sqrt(xb**2+yb**2)
        
        # create new space to interpolate data
        nptNew = self.nptsNew
        nptOld = xb.shape[0]
        x = np.linspace(mesh.xmin,mesh.xmax,num=nptNew)
        y = np.linspace(mesh.ymin,mesh.ymax,num=nptNew)
        XI,YI = np.meshgrid(x,y)
        
        # RBF
        #rbf = interpolate.Rbf(xb, yb, zb, function='inverse',epsilon=2,smooth=0.01)
        #ZI = rbf(XI, YI)
        
        # interp 2d
        #f = interpolate.interp2d(xb, yb, zb, kind='cubic')
        #ZI = f(x,y)
        #ZI = ZI.reshape(100,100)
        
        # ND linear interpolation
        # convert points to 2D array
        xb=xb.reshape(nptOld,1) # need this reshape for hstack to work
        yb=yb.reshape(nptOld,1)
        pt = np.hstack((xb,yb))
        zb=zb.reshape(nptOld,1)
        
        interp = interpolate.LinearNDInterpolator(pt,zb)
        
        XI = XI.flatten()
        YI = YI.flatten()
        XI = XI.reshape(XI.shape[0],1)
        YI = YI.reshape(YI.shape[0],1)
        pt = np.hstack((XI,YI))
        ZI = interp(pt)
        ZI = ZI.reshape(nptNew,nptNew)
        YI = YI.reshape(nptNew,nptNew)
        XI = XI.reshape(nptNew,nptNew)
        
        # 3D figure
        #fig = plt.figure()
        #ax = fig.gca(projection='3d')
        #surf = ax.plot_surface(XI, YI, ZI, cmap=cm.coolwarm,
        #               linewidth=0, antialiased=False)
        
        # contours
        plt.figure()
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.grid(True)
        levels = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
        #levels = np.linspace(0,1.0,num=25)
        cp = plt.contourf(XI,YI,ZI,cmap='hsv',levels=levels)
        #cp = plt.contour(XI,YI,ZI,cmap='hsv',levels=levels)
        plt.clabel(cp,inline=True,fontsize=8)
        #plt.colorbar(cp)
        plt.title('Contour plot')
        
    # march
    def march(self,Nstep,mesh):
        """
        This is the main subroutine that marches the solution in time
        """
        istep  = 0
        deltaT = 1.0E-3
        resplot= 20
        plt.figure(resplot)
        plt.xlabel('Iter')
        plt.ylabel('L2 norm')
        plt.title('Residue plot')
        while istep < Nstep:
            
            self.cf  = np.zeros(self.nn_int)
            self.d2f = np.zeros(self.nn_int)
            self.d4f = np.zeros(self.nn_int)
            R        = np.zeros(self.nn_int)
            
            # calculate switch; undivided difference
            self.sw = np.zeros(self.nn_int)
            self.du = np.zeros(self.nn_int)
            self.__switch(mesh)
            
            self.__calcEulerFlux(mesh)
            self.__calcDissipation2(mesh)
            
            # copy U
            Ucopy = np.zeros(self.nn_int)
            i     = 0
            for u in self.U:
                if u.aux:
                    continue
                Ucopy[i]=u.value
                i       = i + 1
            
            # time step
            res = 0.0
            for irk in range(3):
                R  = -self.cf + self.d2f - self.d4f
                i  = 0
                if irk == 0:
                    R0 = R.copy()
                    for u in self.U:
                         u.value = Ucopy[i] + deltaT*R[i]/u.area
                         i       = i + 1
                if irk > 0:
                    for u in self.U:
                        u.value = Ucopy[i] + 0.5*deltaT*(R[i]+R0[i])/u.area
                        i       = i + 1
            
            i = 0
            for u in self.U:
                res     = res + (u.value-Ucopy[i])**2
                i       = i + 1
            
            
            #self.__BC()
        
            # increment iteration count
            istep = istep + 1
            self.Ctime = istep*deltaT
            
            # print residue
            res = np.sqrt(res / self.nn_int)
            #print("L2norm of delta =",res," after iter =",istep)
            plt.plot(istep,res,'b.')
            plt.pause(0.05)
            
            # plot line
            #if istep%250 == 0 : self.plot_line_graph(mesh,False)
        
    # calculate switch
    def __switch(self,mesh):
        # loop over all interior cells
        for j in range(self.nn_int):
            num = 0.0
            den = 0.0
            for ie in mesh.c2e[j]:
                e = mesh.edges[ie]
                j1 = e.t1
                j2 = e.t2
                if j2 < 0: j2 = -j2
                jn  = j2
                if jn == j: jn = j1
                num = num + np.abs(self.U[jn].value) - np.abs(self.U[j].value)
                den = den + np.abs(self.U[jn].value) + np.abs(self.U[j].value)
                self.du[j] = self.du[j] + self.U[jn].value - self.U[j].value
            
            self.sw[j] = abs(num)/(den+1.0e-6)
            
        
    # end __switch
            
            
    # calculate Euler fluxes
    def __calcEulerFlux(self,mesh):
        Vvec = np.array([self.A,self.B])
        
        # loop over all edges and calculate flux
        for e in mesh.edges:
            j1 = e.t1
            j2 = e.t2
            if j2 < 0:
                j2 = -j2
            ue = 0.5*(self.U[j1].value + self.U[j2].value)
            vc = e.length*np.inner(e.normal,Vvec)
            
            # fluxes; +ve going outward; -ve going inward
            self.cf[j1] = self.cf[j1] + vc*ue
            if e.t2 > 0:
                self.cf[j2] = self.cf[j2] - vc*ue 
    
    # calculate artificial dissipation
    def __calcDissipation2(self,mesh):
        Vvec = np.array([self.A,self.B])
        
        # loop over all edges and calculate flux
        for e in mesh.edges:
            j1 = e.t1
            j2 = e.t2
            if j2 < 0:
                j2 = -j2
            ue = self.U[j2].value - self.U[j1].value
            vc = np.abs(e.length*np.inner(e.normal,Vvec))
            sw = self.k2*max(self.sw[j1],self.sw[j2])
            
            # fluxes; +ve going outward; -ve going inward
            self.d2f[j1] = self.d2f[j1] + sw*vc*ue
            if e.t2 > 0:
                self.d2f[j2] = self.d2f[j2] - sw*vc*ue
            
    
    # 4th order dissipation
    def __calcDissipation4(self,mesh):
        Vvec = np.array([self.A,self.B])
        
        # loop over all edges
        for e in mesh.edges:
            j1 = e.t1
            j2 = e.t2
            if j2 < 0: j2 = - j2
            due = self.du[j2] - self.du[j1]
            vc  = np.abs(e.length*np.inner(e.normal,Vvec))
            sw  = max(0,self.k4-self.k2*self.sw[j1])
            
            # get flux
            self.d4f[j1] = self.d4f[j1] + sw*vc*due
            if e.t2 > 0:
                self.d4f[j2] = self.d4f[j2] - sw*vc*due
            
    # end __calcDissipation4    
    
    
        
    
            
            
            
            
        
