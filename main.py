# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 17:50:36 2020

@author: Shalabh saxena
"""
from my_mesh import mesh
from my_solver import solver

mesh_tri = mesh()
mesh_tri.mesh_domain()
#mesh_tri.show_mesh()

solution = solver(mesh_tri)
#solution.plot_contour(mesh_tri)
linePlot = 5
#plt.figure(linePlot)
solution.plot_line_graph(mesh_tri,linePlot)

# march the solution now
nmax = 2000
solution.march(nmax,mesh_tri)
#solution.plot_contour(mesh_tri)
solution.flagWrite = True
solution.plot_line_graph(mesh_tri,linePlot)
solution.flagWrite = False
