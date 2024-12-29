#from basix.ufl_wrapper import create_vector_element
from basix.ufl import *
from ufl import (Coefficient, Constant, FunctionSpace, Mesh,
                 TestFunction, TrialFunction, ds, dx, grad, inner, triangle, curl, Measure)
'''
shape = "triangle"


degT = 2
degQ = degT

qu_1 = FiniteElement(family = "Quadrature",
                             cell = shape,
                             degree = degQ,
                             quad_scheme="default")

cg_2 = FiniteElement("CG", cell = shape, degree = degT)

qu_2 = VectorElement(qu_1,2)

vt = TestFunction(qu_2)

uq = TrialFunction(cg_2)

c = Constant(shape)

dx = Measure("dx", metadata={"quadrature_degree": degQ})

#a = inner(c*curl(uq), vt) * dx()


'''
shape = "triangle"

dim = 2
degT = 2
degQ = degT

qu_1 = quadrature_element(shape, (dim,), "default", degQ)
cg_2 =  element("CG", shape, degree = degT)

coord_element = element("Lagrange", shape, 1, shape=(2,))
mesh = Mesh(coord_element)

V = FunctionSpace(mesh, qu_1)
solSpace = FunctionSpace(mesh, cg_2)

vt = TestFunction(V)

uq = TrialFunction(solSpace)

c = Constant(mesh)

dx = Measure("dx", metadata={"quadrature_degree": degQ})

a = inner(c*curl(uq), vt) * dx()

