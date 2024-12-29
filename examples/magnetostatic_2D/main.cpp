
//#define ENABLE_DPC_DEBUGGING_FUNCTIONS
#include <magnetics_toolbox/maxwell_solvers.h>
#include <petscsys.h>

//#define ENABLE_DETAILED_DEBUG





int main(int argc, char* argv[]){

    PetscInitialize(&argc, &argv, nullptr, nullptr);    
    {   
        mag_tools::magnetic_field_2D_Az(argc,argv);
    } 

    PetscOptionsSetValue(PETSC_NULLPTR,"-options_left", "no"); 
    PetscFinalize();
    return 0;
    
}