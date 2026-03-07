# magnetix_toolbox
Some calculations of magnetostatic fields based on FEniCSx

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
-    calculation of magnetic forces on bodies using the virtual work method
-    an example calculation of the TEAM Problem 32 including the necessary parameters for the hysteron group model (see Comupumag website for a [general description](https://www.compumag.org/jsite/images/stories/TEAM/problem32.pdf))

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

-   fenicsx 0.10 with
    -   Adios2 

### optional

- dolfinx_mpc 0.10.0 


### python libraries

- matplotlib
- pandas

### optional python libraries
- gmsh 
- meshio

### recommended
- Paraview

# Installation

## Prerequisites

Ensure that FEniCSx libraries version 0.10.0 and their dependencies are available in your environment

## Quick Build (Recommended)

The easiest way to build everything is using the provided build script:

```bash
./build.sh
```

Or using the Python variant:

```bash
python3 build.py
```

Both scripts will:
- Configure and build all components (hysteresis models, fenicsx-magnetics-toolbox, and examples)
- Install everything to the `install/` directory
- Display environment setup instructions

### Build Script Options

```bash
./build.sh -h                 # Show help
./build.sh -j 4               # Use 4 parallel jobs
./build.sh -c                 # Clean build (remove build artifacts)
./build.sh -v                 # Verbose output
```

The same options work with `python3 build.py`.

## Environment Setup

After building, the scripts will display the required environment variables. To make them permanent, add to your shell configuration (`~/.bashrc`, `~/.zshrc`, etc.):

```bash
export CMAKE_PREFIX_PATH="${MAGNETIX_TOOLBOX}/install/dpc_hysteresis-0.1:${MAGNETIX_TOOLBOX}/install/fenicsx_magnetics_toolbox-0.10:$CMAKE_PREFIX_PATH"
export PYTHONPATH="${MAGNETIX_TOOLBOX}/install/fenicsx_magnetics_toolbox-0.10/python:$PYTHONPATH"
```

where `${MAGNETIX_TOOLBOX}` is the path to this repository.

## Manual Build (Alternative)

If you prefer to build manually:

1. **Build hysteresis models**

        cd hysteresis_model
        mkdir build && cd build
        cmake -DCMAKE_INSTALL_PREFIX=../install ..
        make -j 4
        make install

2. **Build fenicsx-magnetics-toolbox**

        cd fenicsx_tools/library
        mkdir build && cd build
        cmake -DCMAKE_INSTALL_PREFIX=../../install -DCMAKE_PREFIX_PATH=../../install/dpc_hysteresis-0.1 ..
        make -j 4
        make install

3. **Build examples**

        cd examples/magnetostatic_2D
        mkdir build && cd build
        cmake -DCMAKE_PREFIX_PATH=../../install/dpc_hysteresis-0.1:../../install/fenicsx_magnetics_toolbox-0.10 ..
        make -j 4

## Notes

- Not having *dolfinx_mpc* installed will result in a warning from cmake which can be ignored.


## TEAM Problem 32

1.  **(Optional)** Generate custom input XML for the example

        cd examples/magnetostatic_2D/TEAM_Problem_32
        python3 TeamProblem32_setup.py

    **Note:** This requires `gmsh` to be installed. Default input scenarios are already prepared and shipped with the repository in the `input/` directory.

2.  Run the example (if possible in parallel)

        cd examples/magnetostatic_2D/TEAM_Problem_32
        mpirun -np 4 ../magnetostatic_2D_exec --scen TeamProblem32_case3_default.xml

3.  All results are found in a time-stamped results directory. 
    If using default inputs: *TEAM_Problem_32_case3_default_<DATE>/results/*
    
    The fields can be viewed e.g. with Paraview. Evaluate flux density at the query points using (remember that **install/fenicsx_magnetics_toolbox-0.10/python** must be on your **PYTHONPATH**)

        python3 evalFieldQuery.py TEAM_Problem_32_case3_default_<DATE>/results/results.xml

# Third party contributions
## tinyXML2
A version of the source code of [tinyXML2](https://github.com/leethomason/tinyxml2) written by Lee Thomason is used and distribued with the hysteresis model.

## msh_to_xdmf
The content of the file *msh_to_xdmf.py* was written by Connor D. Pierce ([see post](https://fenicsproject.discourse.group/t/using-a-simple-3d-mesh-from-gmsh/9639/9)).
