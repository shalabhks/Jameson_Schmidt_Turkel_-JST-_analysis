# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 10:57:25 2020

@author: Shalabh saxena
symbolic manipulator for RK3 on equilateral triangle with only CD expression
checking
"""
import sympy as sym
from sympy.vector import Vector
from sympy.vector import CoordSys3D
sym.init_printing(use_latex=True)

N = CoordSys3D('N')

# symbols
kx,ky      = sym.symbols('kx ky',real=True)
R          = sym.symbols('R',real=True,positive=True)
phix, phiy = sym.symbols('phi_x phi_y',real=True)
A,B        = sym.symbols('A B',real=True)
sga,sgb    = sym.symbols('sigma_a sigma_b',real=True)
U0         = sym.symbols('U_0',real=True)
sg         = sym.symbols('sigma',real=True,positive=True)


# 3 vertices of triangle
v1 = (-sym.sqrt(3)/2)*R*N.i + sym.Rational(-1,2)*R*N.j
v2 = (sym.sqrt(3)/2)*R*N.i + sym.Rational(-1,2)*R*N.j
v3 = R*N.j

# 3 edges
ds1 = v3 - v2
ds2 = v1 - v3
ds3 = v2 - v1

# centers of 3 adjacent triangles
c1 = sym.cos(sym.pi/6)*R*N.i + sym.sin(sym.pi/6)*R*N.j
c2 = sym.cos(sym.pi/6 + 2*sym.pi/3)*R*N.i + sym.sin(sym.pi/6 + 2*sym.pi/3)*R*N.j
c3 = sym.cos(sym.pi/6 + 4*sym.pi/3)*R*N.i + sym.sin(sym.pi/6 + 4*sym.pi/3)*R*N.j

# k vector
k = kx*N.i + ky*N.j

# Area of base triangle
Area = (3*sym.sqrt(3)/4)*R**2

# Scalar at edge mid points with a factor of 2/U0; non-dimentional
Um1_cap = 1 + sym.exp(sym.I*k.dot(c1))
Um2_cap = 1 + sym.exp(sym.I*k.dot(c2))
Um3_cap = 1 + sym.exp(sym.I*k.dot(c3))

# substiture phix for R*kx and phiy for R*ky
Um1_cap = Um1_cap.subs(R*kx,phix)
Um1_cap = Um1_cap.subs(R*ky,phiy)
Um2_cap = Um2_cap.subs(R*kx,phix)
Um2_cap = Um2_cap.subs(R*ky,phiy)
Um3_cap = Um3_cap.subs(R*kx,phix)
Um3_cap = Um3_cap.subs(R*ky,phiy)

# convert them to sin and cos
Um1_cap = Um1_cap.rewrite(sym.sin)
Um2_cap = Um2_cap.rewrite(sym.sin)
Um3_cap = Um3_cap.rewrite(sym.sin)

# dimentional version of Scalar
Um1 = (U0/2)*Um1_cap
Um2 = (U0/2)*Um2_cap
Um3 = (U0/2)*Um3_cap

# Transport velocity vector
V = A*N.i + B*N.j

#Area vectors
S1 = sym.sqrt(3)*c1
S2 = sym.sqrt(3)*c2
S3 = sym.sqrt(3)*c3

# Q(u)is
ds1_x = ds1.dot(N.i)
ds1_y = ds1.dot(N.j)
ds2_x = ds2.dot(N.i)
ds2_y = ds2.dot(N.j)
ds3_x = ds3.dot(N.i)
ds3_y = ds3.dot(N.j)
Qu    = A*Um1*ds1_y + A*Um2*ds2_y + A*Um3*ds3_y
Qu    = Qu - B*Um1*ds1_x - B*Um2*ds2_x - B*Um3*ds3_x
Qu    = Qu/Area
Qu_re = sym.re(Qu)
Qu_im = sym.im(Qu)
Qu_re = Qu_re.simplify()
Qu_im = Qu_im.simplify()
Qu    = Qu_re + sym.I*Qu_im

# Qu from surface integral method
Qu2 = V.dot(S1)*Um1 + V.dot(S2)*Um2 + V.dot(S3)*Um3
Qu2 = Qu2/Area
Qu2_re = sym.re(Qu2)
Qu2_im = sym.im(Qu2)
Qu2_re = Qu2_re.simplify()
Qu2_im = Qu2_im.simplify()
Qu2    = Qu2_re + sym.I*Qu2_im

# The expressions for Qu and Qu2 comes different but on doing (Qu-Qu2).simplify
# we get 0

# we move forward with qu2 because that is how we will make dissipation
# expressions

# substitute simplification assumptions
# 1. replace A/R and B/R by sigma; Note that there should be a dt but that will
# come when we will do alpha*dt
Qu2_re = sym.re(Qu2)
Qu2_im = sym.im(Qu2)
Qu2_re = Qu2_re.subs(A,sg)
Qu2_re = Qu2_re.subs(B,sg)
Qu2_re = Qu2_re*R
Qu2_re = Qu2_re.expand(trig=True)
Qu2_im = Qu2_im.subs(A,sg)
Qu2_im = Qu2_im.subs(B,sg)
Qu2_im = Qu2_im*R
Qu2_im = Qu2_im.expand(trig=True)
Qu2_re = Qu2_re.simplify()
Qu2_im = Qu2_im.simplify()
Qu2    = Qu2_re + sym.I*Qu2_im

"""
*************************2nd Order Diffusion Term******************************
"""

# 2nd order diffusion term D2(U)
dUm1 = U0*(sym.exp(sym.I*k.dot(c1))-1)
dUm2 = U0*(sym.exp(sym.I*k.dot(c2))-1)
dUm3 = U0*(sym.exp(sym.I*k.dot(c3))-1)

# substitute phi_x for kx*R and phi_y for ky*R
dUm1 = dUm1.subs(R*kx,phix)
dUm1 = dUm1.subs(R*ky,phiy)
dUm2 = dUm2.subs(R*kx,phix)
dUm2 = dUm2.subs(R*ky,phiy)
dUm3 = dUm3.subs(R*kx,phix)
dUm3 = dUm3.subs(R*ky,phiy)

# convert them to sin and cos
dUm1 = dUm1.rewrite(sym.sin)
dUm2 = dUm2.rewrite(sym.sin)
dUm3 = dUm3.rewrite(sym.sin)

Du2  = sym.Abs(V.dot(S1))*dUm1 + sym.Abs(V.dot(S2))*dUm2 + sym.Abs(V.dot(S3))*dUm3
Du2  = Du2/Area
Du2_re = sym.re(Du2)
Du2_im = sym.im(Du2)
# sub Adt/R and Bdt/R by sg
Du2_re = Du2_re.subs(A,sg)
Du2_re = Du2_re.subs(B,sg)
Du2_re = Du2_re*R
Du2_re = Du2_re.expand(trig=True)
Du2_im = Du2_im.subs(A,sg)
Du2_im = Du2_im.subs(B,sg)
Du2_im = Du2_im*R
Du2_im = Du2_im.expand(trig=True)
Du2_re = Du2_re.simplify()
Du2_im = Du2_im.simplify()
Du2    = Du2_re + sym.I*Du2_im
#display(Du2_re)
#display(Du2_im)

"""
***************************4th Order Dissipation*******************************
"""
delSQ_U0 = dUm1 + dUm2 + dUm3
delSQ_U0 = delSQ_U0.expand(trig=True)
delSQ_U0 = delSQ_U0.simplify()
delSQ_U0_re = sym.re(delSQ_U0)
delSQ_U0_im = sym.im(delSQ_U0)
delSQ_U0 = delSQ_U0_re + sym.I*delSQ_U0_im
#display(delSQ_U0)

# del square U1
delSQ_U1 = (delSQ_U0)*sym.exp(sym.I*k.dot(c1))
delSQ_U1 = delSQ_U1.simplify()
delSQ_U1_re = sym.re(delSQ_U1)
delSQ_U1_im = sym.im(delSQ_U1)

delSQ_U1_re = delSQ_U1_re.subs(kx*R,phix)
delSQ_U1_re = delSQ_U1_re.subs(ky*R,phiy)
delSQ_U1_re = delSQ_U1_re.simplify()
delSQ_U1_re = delSQ_U1_re.expand(trig=True)
delSQ_U1_re = delSQ_U1_re.collect(U0)

delSQ_U1_im = delSQ_U1_im.subs(kx*R,phix)
delSQ_U1_im = delSQ_U1_im.subs(ky*R,phiy)
delSQ_U1_im = delSQ_U1_im.simplify()
delSQ_U1_im = delSQ_U1_im.expand(trig=True)
delSQ_U1_im = delSQ_U1_im.collect(U0)

delSQ_U1 = delSQ_U1_re + sym.I*delSQ_U1_im

# del square U2
delSQ_U2 = (delSQ_U0)*sym.exp(sym.I*k.dot(c2))
delSQ_U2 = delSQ_U2.simplify()
delSQ_U2_re = sym.re(delSQ_U2)
delSQ_U2_im = sym.im(delSQ_U2)

delSQ_U2_re = delSQ_U2_re.subs(kx*R,phix)
delSQ_U2_re = delSQ_U2_re.subs(ky*R,phiy)
delSQ_U2_re = delSQ_U2_re.simplify()
delSQ_U2_re = delSQ_U2_re.expand(trig=True)
delSQ_U2_re = delSQ_U2_re.collect(U0)

delSQ_U2_im = delSQ_U2_im.subs(kx*R,phix)
delSQ_U2_im = delSQ_U2_im.subs(ky*R,phiy)
delSQ_U2_im = delSQ_U2_im.simplify()
delSQ_U2_im = delSQ_U2_im.expand(trig=True)
delSQ_U2_im = delSQ_U2_im.collect(U0)

delSQ_U2 = delSQ_U2_re + sym.I*delSQ_U2_im

# del square U3
delSQ_U3 = (delSQ_U0)*sym.exp(sym.I*k.dot(c3))
delSQ_U3 = delSQ_U3.simplify()
delSQ_U3_re = sym.re(delSQ_U3)
delSQ_U3_im = sym.im(delSQ_U3)

delSQ_U3_re = delSQ_U3_re.subs(kx*R,phix)
delSQ_U3_re = delSQ_U3_re.subs(ky*R,phiy)
delSQ_U3_re = delSQ_U3_re.simplify()
delSQ_U3_re = delSQ_U3_re.expand(trig=True)
delSQ_U3_re = delSQ_U3_re.collect(U0)

delSQ_U3_im = delSQ_U3_im.subs(kx*R,phix)
delSQ_U3_im = delSQ_U3_im.subs(ky*R,phiy)
delSQ_U3_im = delSQ_U3_im.simplify()
delSQ_U3_im = delSQ_U3_im.expand(trig=True)
delSQ_U3_im = delSQ_U3_im.collect(U0)

delSQ_U3 = delSQ_U3_re + sym.I*delSQ_U3_im

rex = delSQ_U1_re - delSQ_U0_re
imx = delSQ_U1_im - delSQ_U0_im
rex = rex.simplify()
imx = imx.simplify()
d_delSQ_U1 = rex + sym.I*imx

rex = delSQ_U2_re - delSQ_U0_re
imx = delSQ_U2_im - delSQ_U0_im
rex = rex.simplify()
imx = imx.simplify()
d_delSQ_U2 = rex + sym.I*imx

rex = delSQ_U3_re - delSQ_U0_re
imx = delSQ_U3_im - delSQ_U0_im
rex = rex.simplify()
imx = imx.simplify()
d_delSQ_U3 = rex + sym.I*imx

Du4  = sym.Abs(V.dot(S1))*d_delSQ_U1 + sym.Abs(V.dot(S2))*d_delSQ_U2 +        \
sym.Abs(V.dot(S3))*d_delSQ_U3
Du4  = Du4/Area
Du4_re = sym.re(Du4)
Du4_im = sym.im(Du4)
# sub Adt/R and Bdt/R by sg
Du4_re = Du4_re.subs(A,sg)
Du4_re = Du4_re.subs(B,sg)
Du4_re = Du4_re*R
Du4_re = Du4_re.expand(trig=True)
Du4_im = Du4_im.subs(A,sg)
Du4_im = Du4_im.subs(B,sg)
Du4_im = Du4_im*R
Du4_im = Du4_im.expand(trig=True)
Du4_re = Du4_re.simplify()
Du4_im = Du4_im.simplify()
Du4_im = Du4_im*3/(4*U0*sg)
Du4_im = Du4_im.expand(trig=True)
Du4_im = Du4_im * (4*U0*sg)/3
Du4    = Du4_re + sym.I*Du4_im

"""
# ut + alpha*u = 0
Qu    = Qu.subs(A,sga)
Qu    = Qu.subs(B,sgb)
alpha_dt = 2*Qu/(3*sym.sqrt(3)*R) # because Qu was already scaled by 2*u0
alpha_dt = alpha_dt.simplify()
"""
