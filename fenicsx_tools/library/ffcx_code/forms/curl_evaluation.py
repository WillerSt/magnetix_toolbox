# SPDX-FileCopyrightText: 2025 Stephan Willerich
# SPDX-License-Identifier: MIT License

#from basix.ufl_wrapper import create_vector_element
from basix.ufl import *
from ufl import (Coefficient, Constant, FunctionSpace, Mesh,
                 TestFunction, TrialFunction, ds, dx, grad, inner, triangle, curl, Measure)

dim = 2

if dim == 2:
    shape = "triangle"
elif dim == 3:
    shape = "tetrahedron"


degT = 2
degQ = degT

#### Definition of Elements
qu_1 = quadrature_element(shape, (dim,), "default", degQ)

quScaElem = quadrature_element(shape, (), "default", degQ)

quVecElem = quadrature_element(shape, (dim,), "default", degQ)
#quVecElem = VectorElement(quScaElem)

#qu_2 = (qu_1,2)

if dim == 2:
    cg_2 =  element("CG", shape, degree = degT)
elif dim == 3:
    cg_2 = element("N1curl", shape, degree = degT )


coord_element = element("Lagrange", shape, 1, shape=(2,))
mesh = Mesh(coord_element)

V = FunctionSpace(mesh, quVecElem)
solSpace = FunctionSpace(mesh, cg_2)

#### Definition of Functions for the calculations on quad points
v = TestFunction(V)
u = TrialFunction(V)

sol =  Coefficient(solSpace)

c = Constant(mesh)

dx = Measure("dx", metadata={"quadrature_degree": degQ})

a = inner(u, v) * dx()
L = inner(c * curl(sol), v) * dx()
'''

shape = "triangle"

degT = 2
degQ = degT

#### Definition of Elements
qu_1 = FiniteElement(family = "Quadrature",
                             cell = shape,
                             degree = degQ,
                             quad_scheme="default")

qu_2 = VectorElement(qu_1,2)

cg_2 =  FiniteElement("CG", cell = shape, degree = degT)


#### Definition of Functions for the calculations on quad points
v = TestFunction(qu_2)
u = TrialFunction(qu_2)

sol =  Coefficient(cg_2)

c = Constant(shape)

dx = Measure("dx", metadata={"quadrature_degree": degQ})

a = inner(u, v) * dx()
L = inner(c * curl(sol), v) * dx()
'''