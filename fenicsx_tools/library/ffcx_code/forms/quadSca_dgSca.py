#from basix.ufl_wrapper import create_vector_element
from basix.ufl import *
from ufl import (Coefficient, Constant, FunctionSpace, Mesh,
                 TestFunction, TrialFunction, ds, dx, grad, inner, triangle, curl, Measure)

shape = "triangle"

dim = 1
degT = 1
degQ = 2

#### Definition of Elements
quElem = quadrature_element(shape, (), "default", degQ)

dgElem =  element("DG", shape, degT)

coord_element = element("Lagrange", shape, 1, shape=(2,))
mesh = Mesh(coord_element)

dg = FunctionSpace(mesh, dgElem)
qu= FunctionSpace(mesh, quElem)


#### Definition of Functions for the calculations on quad points
v = TestFunction(dg)
u = TrialFunction(dg)

sol =  Coefficient(qu)

dx = Measure("dx", metadata={"quadrature_degree": degQ})

a = inner(u, v) * dx()
L = inner(sol, v) * dx()
