#from basix.ufl_wrapper import create_vector_element
from basix.ufl import *
from ufl import (Coefficient, Constant, FunctionSpace, Mesh,
                 TestFunction, TrialFunction, ds, dx, grad, inner, triangle, curl, Measure)

shape = "triangle"

dim = 2
degT = 1
degQ = 2

#### Definition of Elements
quElem = quadrature_element(shape, (dim,), "default", degQ)

quMatElem = quadrature_element(shape, (dim,dim), "default", degQ)

dgElem =  element("DG", shape, degT, shape=(2,))

coord_element = element("Lagrange", shape, 1, shape=(2,))
mesh = Mesh(coord_element)

dg_2 = FunctionSpace(mesh, dgElem)
qu_2= FunctionSpace(mesh, quElem)

qu_mat = FunctionSpace(mesh, quMatElem)




#### Definition of Functions for the calculations on quad points
v = TestFunction(dg_2)
u = TrialFunction(dg_2)

rotFunc = Coefficient(qu_mat)

sol =  Coefficient(qu_2)

dx = Measure("dx", metadata={"quadrature_degree": degQ})

a = inner(u, v) * dx()
L = inner(rotFunc*sol, v) * dx()
