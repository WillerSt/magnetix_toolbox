# magnetix_toolbox
Some calculations of magnetostatic fields, based on FEniCSx

The code was originally written to test hysteresis models and was afterwards enhanced to be of more general use.
In its current state the library only supports 2D calculations based on the z-component of the magnetic vector potential. 


What is working so far:
-    coupling of various magnetic material models to the fenicsx framework:
      -    hysteresis models (hysteron group model is implemented, described in this [paper](https://doi.org/10.1109/TMAG.2019.2954580))
      -    spline interpolation
      -    atan curve
      -    linear materials
-    field sources:
    -    prescribed current density
    -    linear permanent magnet
-    calculation magnetic forces on bodies using the virtual work method
-    an example calculation of the TEAM Problem 32 including the necessary parameters for the hysteron group model

What could be improved:
-    unpolished code with some stubs and deprecated routines
-    no python api (interfacing with material models has to happen on the c++-level)
-    more examples are on the way 

# Prerequisites
## install requirements

- hysteresis_model:
    -   Eigen3
    -   boost_thread  
    -   boost_system
    -   boost_filesystem 

-   fenicsx 0.9 with
    -   Adios2 

### optional

- dolfinx_mpc 0.9.0 

## example requirements

- gmsh 

### python libraries
- meshio
- matplotlib
- pandas

### recommended
- Paraview

# Installation
## hysteresis models
The library containing the hysteron group model needed to run the *TEAM Problem 32* is part of the repository
Install first using

    cd hysteresis_model
    mkdir build && cd build
    cmake ..
    make install

The created library should now be available in the install *directory* located in the base directory of the repository.

## fenicsx-magnetics-toolbox
Make sure that cmake can find the hysteresis model libraries by appending *install/dpc_hysteresis-0.1* to your *CMAKE_PREFIX_PATH*.

Not having *dolfin_mpc* installed installed will result in a warning from cmake which can be ignored.

Install by running

    cd fenicsx_tools/library
    mkdir build && cd build
    cmake ..
    make install

The created library should now be available in the *install* directory located in the base directory of the repository.


# Running the examples
## Preparation
-    Append *install/dpc_hysteresis-0.1*  and the *install/fenicsx_magnetics_toolbox-0.9* directories (default installation paths) to your *CMAKE_PREFIX_PATH*
-    Append *install/fenicsx_magnetics_toolbox-0.9/python/* to your *PYTHON_PATH*


## TEAM Problem 32

1.  create the executable for the 2D magnetostatic field problem
    
        cd  examples/magnetostatic_2D
        mkdir build &&  cd build 
        cmake ..
        make 
        mv magnetostatic_2D_exec .. && cd ..

2. create input xml for example

        cd TEAM_Problem_32
        python TeamProblem32_setup.py
        cd ..

3. run example (if possible in parallel)

       mpirun -np 4 ./magnetostatic_2D_exec --scen TEAM_Problem_32/TeamProblem32_case3.xml

All results are found in the *TEAM_Problem_32/results/* directory. The fields can be viewed e.g. with Paraview.
An automated evaluation of the results is currently being worked on.

# Third party contributions
## tinyXML2
A version of the source code of [tinyXML2](https://github.com/leethomason/tinyxml2) written by Lee Thomason is used and distribued with the hysteresis model.

## msh_to_xdmf
The content of the file *msh_to_xdmf.py* was written by Connor D. Pierce ([see post](https://fenicsproject.discourse.group/t/using-a-simple-3d-mesh-from-gmsh/9639/9)).
