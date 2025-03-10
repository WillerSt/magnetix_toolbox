# SPDX-FileCopyrightText: 2025 Stephan Willerich
# SPDX-License-Identifier: MIT License

from basix.ufl import *
from ufl import (Coefficient, Constant, FunctionSpace, Mesh,
                 TestFunction, TrialFunction, ds, dx, grad, inner, triangle, curl, Measure)

shape = "triangle"

dim = 2
degT = 1
degQ = 2

#### Definition of Elements
quElem = quadrature_element(shape, (dim,), "default", degQ)

dgElem =  element("DG", shape, degT, shape=(2,))

coord_element = element("Lagrange", shape, 1, shape=(2,))
mesh = Mesh(coord_element)

dg_2 = FunctionSpace(mesh, dgElem)
qu_2= FunctionSpace(mesh, quElem)


#### Definition of Functions for the calculations on quad points
v = TestFunction(dg_2)
u = TrialFunction(dg_2)

sol =  Coefficient(qu_2)

dx = Measure("dx", metadata={"quadrature_degree": degQ})

a = inner(u, v) * dx()
L = inner(sol, v) * dx()
