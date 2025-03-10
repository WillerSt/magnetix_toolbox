# SPDX-FileCopyrightText: 2025 Stephan Willerich
# SPDX-License-Identifier: MIT License

from basix.ufl import *
from ufl import (Coefficient, Constant, FunctionSpace, Mesh, form,
                 TestFunction, TrialFunction, dot, ds, dx, grad, inner, triangle, curl, Measure, replace, derivative, SpatialCoordinate)

shape = "triangle"

degT = 2
degQ = degT
dim = 2

coord_element = element("Lagrange", shape, 1, shape=(dim,))
msh = Mesh(coord_element)

quScaElem = quadrature_element(shape, (), "default", degQ)
qu_1 = FunctionSpace(msh, quScaElem)

quVecElem = quadrature_element(shape, (dim,), "default", degQ)
qu_2 = FunctionSpace(msh, quVecElem)

quMatElem = quadrature_element(shape, (dim,dim), "default", degQ)
qu_jacobean = FunctionSpace(msh, quMatElem)

dg0Elem  =  element("DG", shape, 0)
dg_0 = FunctionSpace(msh, dg0Elem)

cgElem=  element("CG", shape, degT)
cg_2 = FunctionSpace(msh, cgElem)

A = Coefficient(cg_2)
B  = Coefficient(qu_2)
H  = Coefficient(qu_2)
nu = Coefficient(qu_jacobean)
vacInd = Coefficient(dg_0)

Wmag0 = Coefficient(qu_1)

dx= Measure("dx", metadata={"quadrature_degree": degQ})

#Wmag = 1/2*dot(B,H)*dx
#dWmag = dot(H,curl(A))*dx() - dot(nu*B,curl(A))*dx() + 1/2*dot(nu*curl(A),curl(A))*dx()

#Wmag = Wmag0*dx() + dot(H,curl(A))*dx() - dot(nu*B,curl(A))*dx() + 1/2*dot(nu*curl(A),curl(A))*dx() + 1/2*dot(nu*B,B)*dx() - dot(H,B)*dx()

c0 = Constant(msh)
c1 = Constant(msh)
c2 = Constant(msh)
Wmag = c0*Wmag0*dx() + c1*1/2*dot(nu*B,B)*dx() - c1*dot(H,B)*dx() + c2*dot(H,curl(A))*dx() - c2*dot(nu*B,curl(A))*dx()  + vacInd * 1/2*dot(nu*curl(A),curl(A))*dx() 
X = SpatialCoordinate(msh)
P1 = FunctionSpace(msh, coord_element)
v = TestFunction(P1)
L = derivative(Wmag, X, v)







