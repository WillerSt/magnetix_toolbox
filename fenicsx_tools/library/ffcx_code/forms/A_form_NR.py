# SPDX-FileCopyrightText: 2025 Stephan Willerich
# SPDX-License-Identifier: MIT License

#from basix.ufl_wrapper import create_vector_element
from basix.ufl import *
from ufl import (Coefficient, Constant, FunctionSpace, Mesh,
                 TestFunction, TrialFunction, ds, dx, grad, inner, triangle, curl, Measure)

#### Basic Definitions of Shape and degree
# Compiles with:  ffc -O -fquadrature_degree=2 -fsplit -v -l dolfin A_form_NR.ufl

shape = "triangle"

degT = 2
degQ = degT
dim = 2

#### Definition of Elements

quScaElem = quadrature_element(shape, (), "default", degQ)

quVecElem = quadrature_element(shape, (dim,), "default", degQ)
#quVecElem = VectorElement(quScaElem)

quMatElem = quadrature_element(shape, (dim,dim), "default", degQ)
#quMatElem = TensorElement(quScaElem,2)

cgElem=  element("CG", shape, degT)

coord_element = element("Lagrange", shape, 1, shape=(2,))
mesh = Mesh(coord_element)


##### Definition of Functions needed for the solved form

qu_1 = FunctionSpace(mesh, quScaElem)
qu_2 = FunctionSpace(mesh, quVecElem)
qu_jacobean = FunctionSpace(mesh, quMatElem)
cg_2 = FunctionSpace(mesh, cgElem)

M = Coefficient(qu_2)

v = TestFunction(cg_2)
A_delta =  TrialFunction(cg_2)

A_old = Coefficient(cg_2)

nu = Constant(mesh)

J = Coefficient(qu_jacobean)
j = Coefficient(qu_1)
M0 = Coefficient(qu_2)

dx = Measure("dx", metadata={"quadrature_degree": degQ})

a =inner(J * curl(A_delta), curl(v)) * dx()

L = inner(j, v)*dx() + inner(M0, curl(v)) * dx() + inner(M, curl(v)) * dx() - nu * inner(curl(A_old),curl(v)) * dx()

'''
shape = "triangle"

degT = 2
degQ = degT

#### Definition of Elements
cg_2 =  FiniteElement("CG", cell = shape, degree = degT)

qu_1 = FiniteElement(family = "Quadrature",
                             cell = shape,
                             degree = degQ,
                             quad_scheme="default")
qu_2 = VectorElement(qu_1,2)


qu_jacobean = TensorElement(qu_1,2)

##### Definition of Functions needed for the solved form


M = Coefficient(qu_2)

v = TestFunction(cg_2)
A_delta =  TrialFunction(cg_2)

A_old = Coefficient(cg_2)

nu = Constant(shape)

J = Coefficient(qu_jacobean)
j = Coefficient(qu_1)

dx = Measure("dx", metadata={"quadrature_degree": degQ})

a =inner(J * curl(A_delta), curl(v)) * dx()

L = inner(j, v)*dx() + inner(M, curl(v)) * dx() - nu * inner(curl(A_old),curl(v)) * dx()

'''