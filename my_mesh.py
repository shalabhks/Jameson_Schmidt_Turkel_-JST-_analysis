# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 16:31:16 2020

@author: Shalabh saxena
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import spatial

class edge_ds:
    v1    = -1             # vertex 1
    v2    = -1             # vertex 2; v1 < v2
    t1    = -1             # triangle 1
    t2    = -1             # triangle 2
    bndry = -1             # 0 -> interior; 1->boundary
    center = np.array([2]) 
    length = 0.0
    normal = np.array([2])


class mesh:
    points    = np.array([])
    tri       = np.array([])
    c2e       = {}
    nedges    = 0
    nBCedges  = 0
    npoints   = 0
    nelem     = 0
    nelemInt  = 0
    edges     = []
    xmin      = -1.0
    xmax      =  1.0
    ymin      = -1.0
    ymax      =  1.0
    
    # create mesh
    def mesh_domain(self):
        ny = 50
        nx = 50
        
        # create points inside internal domain
        self.points  = np.zeros([nx*ny,2])
        self.npoints = nx*ny
        yd = np.linspace(-1.0,1.0,ny)
        xd = np.linspace(-1.0,1.0,nx)
        k  = 0
        for i in range(nx):
            for j in range(ny):
                xx = xd[i]
                yy = yd[j]
                self.points[k,0] = xx
                self.points[k,1] = yy
                
                k = k + 1
         
        print("\n Total points =",self.npoints)
        
        #self.__plot_points()
        my_pts    = self.points.copy()
        my_tri = spatial.Delaunay(my_pts,qhull_options="QJ")
        my_tri.close
        
        # remove slivers if any
        self.tri = self.__remove_slivers(my_tri.simplices,my_pts)
        
        # fill total number of tri
        self.nelem    = self.tri.shape[0]
        self.nelemInt = self.nelem
        
        # create edge d/s
        self.__create_edges()
        
        # create periodic bc connections
        self.__assign_periodic()
        # fill data structure
        #self.tri  = my_tri.simplices
        
        # smooth mesh
        self.__smoothGraph()
        
    # plot points
    def __plot_points(self):
        fig = plt.figure()
        plt.plot(self.points[:,0],self.points[:,1],'.')
    
    # this subroutine plots the mesh    
    def show_mesh(self):
        plt.figure(figsize=(6,6))
        plt.triplot(self.points[:,0],self.points[:,1],self.tri)
    

    # this subroutine removes slivers from triangles
    def __remove_slivers(self,my_tri,my_pts):
        cnt = 0
        ok  = []
        j   = 0
        for tri in my_tri:
            pts = my_pts[tri]
            v1  = pts[1,:] - pts[0,:]
            v2  = pts[2,:] - pts[0,:]
            cr  = np.cross(v1,v2)
            ar  = 0.5*np.sqrt(np.inner(cr,cr))
            j   = j + 1
            if ar <= 0.0:
                cnt =cnt + 1
                next
            else:
                ok.append(j-1)
        
            
        print("# total elements        =",my_tri.shape[0])
        print("# Sliver elements found =",cnt)
        print("size of ok =",len(ok))
        my_tri2 = my_tri[ok]
        
        return my_tri2
                
    # This subroutine creates edge d/s
    def __create_edges(self):
        my_tri=self.tri.copy()
        
        # empty graph (dictionary); initialize with only keys and empty lists
        graph = dict.fromkeys(range(self.points.shape[0]),[])
        
        """
There is probably a python bug here. If lst is not used then all entries of
graph are appended in each of the below statement and not only with the given 
key.
        """
        for t in my_tri:
            # extend graph for vertex 0
            lst = graph[t[0]].copy()
            try:
                lst.index(t[1])
            except ValueError:
                lst.append(t[1])
                
            try:
                lst.index(t[2])
            except ValueError:
                lst.append(t[2])
            graph[t[0]]=lst
            
            # extend graph for vertex 1
            lst = graph[t[1]].copy()
            try:
                lst.index(t[0])
            except ValueError:
                lst.append(t[0])
                
            try:
                lst.index(t[2])
            except ValueError:
                lst.append(t[2])
            graph[t[1]]=lst
            
            # extend graph for vertex 2
            lst = graph[t[2]].copy()
            try:
                lst.index(t[0])
            except ValueError:
                lst.append(t[0])
                
            try:
                lst.index(t[1])
            except ValueError:
                lst.append(t[1])
            graph[t[2]]=lst
            
        # report graph statistics;# we do not remove orphaned points
        cnt0 = 0
        for i in range(self.points.shape[0]):
            #print("kkkk",graph[i])
            if len(graph[i]) == 0:
                cnt0 = cnt0 + 1
            
        if cnt0 != 0:
            print("\n Numer of orphaned points =",cnt0)
        
        # create edge d/s; v1,v2,length
        for i in range(self.points.shape[0]):
            for vv in graph[i]:
                if vv < i:
                    next
                e        = edge_ds()
                e.v1     = i
                e.v2     = vv
                dv       = self.points[i]-self.points[vv]
                e.length = np.sqrt(np.inner(dv,dv))
                e.center = 0.5*(self.points[i]+self.points[vv])
                self.edges.append(e)
                # remove vertex i from list of vertex vv
                graph[vv].remove(i)
            
        self.nedges = len(self.edges)
        print("#edges =",self.nedges)
        del graph
                
        # Now we create a list of elems associated with each vertex
        elem2points = dict.fromkeys(range(self.points.shape[0]),[])
        i           = 0
        for t in my_tri:
            # to first vertex
            lst = elem2points[t[0]].copy()
            lst.append(i)
            elem2points[t[0]] = lst
            
            # to 2nd vertex
            lst = elem2points[t[1]].copy()
            lst.append(i)
            elem2points[t[1]] = lst
            
            # to 3rd vertex
            lst = elem2points[t[2]].copy()
            lst.append(i)
            elem2points[t[2]] = lst
            
            # increment i
            i = i + 1
        
        # also create a list of edges to vert
        edge2points = dict.fromkeys(range(self.nedges),[])
        i           = 0
        for e in self.edges:
            lst = edge2points[e.v1].copy()
            lst.append(i)
            edge2points[e.v1]=lst
            
            lst = edge2points[e.v2].copy()
            lst.append(i)
            edge2points[e.v2]=lst
            
            i = i + 1
        
        # now loop through all elements and for each edge find list of elems
        # attached to them. One of the common elem will be this elem, find other
        # one. That is e2. In case there is no other common elem then it is
        # boundary edge
        self.c2e = dict.fromkeys(range(self.tri.shape[0]),[])
        nbc_edge = 0
        i        = 0
        for t in my_tri:
            lst = []
            for j in range(3):
                jn    = ((j%3)+1)%3  # Ha Ha Ha
                Ledv1 = set(edge2points[t[j ]].copy()) # as set not list
                Ledv2 = set(edge2points[t[jn]].copy())
                Lelv1 = set(elem2points[t[j ]].copy())
                Lelv2 = set(elem2points[t[jn]].copy())
            
                ie = list(Ledv1.intersection(Ledv2)) # edge id
                ie = ie[0]      # list to scalar
                lst.append(ie)  # append to c2e list
                ce = list(Lelv1.intersection(Lelv2)) # sharing elements
            
                if len(ce)==1:
                    nbc_edge           = nbc_edge + 1
                    self.edges[ie].t1  = ce[0]
                    self.edges[ie].bndry = 1 # boundary
                    # sanity check
                    if ce[0] != i:
                        print("Error identifying attatched element : type 1")
                elif len(ce)==2:
                    self.edges[ie].t1 = min(ce)
                    self.edges[ie].t2 = max(ce)
                    self.edges[ie].bndry = 0 # interior
                    # sanity check
                    if ce[0] != i and ce[1] != i:
                        print("Error identifying attatched element : type 2")
                else:
                    print("Error identifying attatched element : type 3")
                
            self.c2e[i] = lst
            i           = i + 1
        
        self.nBCedges = nbc_edge
        print("# Boundary edges =",self.nBCedges)
              
        # fill in the normals
        for e in self.edges:
            # ce3ll center
            p1 = self.points[self.tri[e.t1]] # 3 points
            c1 = np.array([np.mean(p1[:,0]),np.mean(p1[:,1])])
            
            # normal to edge
            p1 = self.points[e.v1]
            p2 = self.points[e.v2]
            v1 = p2-p1
            n1 = np.array([v1[1],-v1[0]])
            
            # points outward
            if np.inner(n1,(e.center-c1)) < 0:
                n1 = - n1
            
            e.normal = n1/e.length
    
                
    # assign periodic connectivity to boundary edges
    def __assign_periodic(self):
        """
Periodic connections are made here. Boundary edges till this place have j2=-1
Here they are set to -ve of element number they are to take their value from.
        """
        xperiod = self.xmax - self.xmin
        yperiod = self.ymax - self.ymin
        
        xmin_e = []
        xmax_e = []
        ymin_e = []
        ymax_e = []
        dx     = xperiod/self.nelem
        dy     = yperiod/self.nelem
        i      = 0
        for e in self.edges:
            if e.bndry < 1:
                i = i + 1
                continue
            if e.center[0] < (self.xmin + dx):
                xmin_e.append(i)
            elif e.center[0] > (self.xmax - dx):
                xmax_e.append(i)
            elif e.center[1] < (self.ymin + dy):
                ymin_e.append(i)
            elif e.center[1] > (self.ymax - dy):
                ymax_e.append(i)
            else:
                print("Error >> couldn't assign boundary edge to side!")
            
            i = i + 1
        
        # set the periodic connectivity along x axis
        #self.show_mesh()
        yy1 = np.array([])
        for i in xmin_e:
           yy1 = np.append(yy1,self.edges[i].center[1])
        yy2 = np.array([])
        for i in xmax_e:
            yy2 = np.append(yy2,self.edges[i].center[1])
        arg1 = np.argsort(yy1)
        arg2 = np.argsort(yy2)
 
        # assumption 1-1 connectivity
        for i in range(len(yy1)):
            ie1 = xmin_e[arg1[i]]
            ie2 = xmax_e[arg2[i]]
            if self.edges[ie1].t2 != -1:
                print("Something wrong! #1")
            if self.edges[ie2].t2 != -1:
                print("Something wrong! #2")
            self.edges[ie1].t2 = -self.edges[ie2].t1
            self.edges[ie2].t2 = -self.edges[ie1].t1
        
        # set the periodic connectivity along y axis
        xx1 = np.array([])
        for i in ymin_e:
           xx1 = np.append(xx1,self.edges[i].center[0])
        xx2 = np.array([])
        for i in ymax_e:
            xx2 = np.append(xx2,self.edges[i].center[0])
        arg1 = np.argsort(xx1)
        arg2 = np.argsort(xx2)
 
        # assumption 1-1 connectivity
        for i in range(len(xx1)):
            ie1 = ymin_e[arg1[i]]
            ie2 = ymax_e[arg2[i]]
            if self.edges[ie1].t2 != -1:
                print("Something wrong! #1")
            if self.edges[ie2].t2 != -1:
                print("Something wrong! #2")
            self.edges[ie1].t2 = -self.edges[ie2].t1
            self.edges[ie2].t2 = -self.edges[ie1].t1
        
        # check
        """
        i = 0
        for e in self.edges:
            if e.t2 >= -1:
                continue
            self.__plotTri(e.t1)
            self.__plotTri(-e.t2)
            i = i + 1
         """
         
    def __smoothGraph(self):
        """
        all boundary nodes are fixed. All internal nodes move
        """
        print("Smoothing mesh---")
        onBoundary     = 1
        nextToBoundary = 2
        Interior       = 0
        
        graph = dict.fromkeys(range(self.points.shape[0]),[])
        ptb   = np.zeros(self.npoints)  # mark boundary points
        
        # build graph and mark boundary nodes
        for e in self.edges:
            lst = graph[e.v1].copy()
            lst.append(e.v2)
            graph[e.v1] = lst
            lst = graph[e.v2].copy()
            lst.append(e.v1)
            graph[e.v2] = lst
            if e.bndry == onBoundary: # boundary edge
                ptb[e.v1] = onBoundary
                ptb[e.v2] = onBoundary
        
        # mark nodes next to boundary
        for i in range(self.npoints):
            if ptb[i]==onBoundary:
                lst=graph[i].copy()
                for j in lst:
                    if ptb[j]==Interior:
                        ptb[j]=nextToBoundary  # next to boundary
        
        # loop to smooth
        res = 1.0
        itr = 0
        while res > 1.0e-6:
            itr = itr + 1
            
            i   = 0
            res = 0.0
            for pt in self.points:
                if ptb[i]==onBoundary:
                    i = i + 1
                    continue
                if ptb[i]==nextToBoundary:
                    coef = 0.25
                if ptb[i]==Interior:
                    coef = 0.5
                
                # find new point by mean
                lst = graph[i].copy()
                ptn = np.zeros((1,2))
                npt = len(lst)
                for j in lst:
                    ptn = ptn + self.points[j]
                ptn = ptn / npt
                
                # update point and res
                delt           = coef*(ptn-pt)
                res            = res + np.inner(delt,delt)
                ptn            = pt + delt
                self.points[i] = ptn
                i              = i + 1
                
            # print residue to screen
            res = np.sqrt(res/self.npoints)
            #print("Residue = ",res," # iter = ",itr)
        
        print("Smoothing completed") 
            
            
    def __plotTri(self,indx):
        t  = self.tri[indx]
        xx = self.points[t,0]
        yy = self.points[t,1]
        plt.plot(xx,yy,'r.')
                    
            
            
            
            
            
            
            
        
        
            
        
                
        
                
                
                
        
