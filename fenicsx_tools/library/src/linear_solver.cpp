#include <magnetics_toolbox/linear_solver.h>
#include <dolfinx/fem/petsc.h>

namespace mag_tools{
    template<typename T> void linear_solver<T>::assemble_matrix(){
        MatZeroEntries(A.mat());
        dolfinx::fem::assemble_matrix(la::petsc::Matrix::set_block_fn(A.mat(), ADD_VALUES),
                        *a, bc);
        MatAssemblyBegin(A.mat(), MAT_FLUSH_ASSEMBLY);
        MatAssemblyEnd(A.mat(), MAT_FLUSH_ASSEMBLY);
            
        dolfinx::fem::set_diagonal<T>(dolfinx::la::petsc::Matrix::set_fn(A.mat(), INSERT_VALUES), 
                        *(a->function_spaces()[0]),
                        bc);
        MatAssemblyBegin(A.mat(), MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(A.mat(), MAT_FINAL_ASSEMBLY);

    }

    template<typename T> void linear_solver<T>::assemble_rhs(){
        this->assemble_rhs_internally(this->L);
    }

    template<typename T> double linear_solver<T>::calc_orig_rhs_norm(const std::shared_ptr<const varForm<T>>& zeroL){
        // assemble_rhs() has to be called externally to get back to old state
        assemble_rhs_internally(zeroL);
        return get_rhs_norm();         
    }

    template<typename T> double linear_solver<T>::get_rhs_norm() const{
        return dolfinx::la::norm(yForm);
        /*
        double norm;
        y.set_from_options();
        VecNorm(y.vec(), NORM_2, &norm);
        

        return double(norm);
        */
    }

    template<typename T> 
    void linear_solver<T>::assemble_rhs_internally(const std::shared_ptr<const varForm<T>>& formIn){
        yForm.set(0.0);
        dolfinx::fem::assemble_vector(yForm.mutable_array(), *formIn);
        dolfinx::fem::apply_lifting<T,U<T>>(yForm.mutable_array(), {this->a}, {this->bc}, {}, T(1.0));
        yForm.scatter_rev(std::plus<T>());

        
        //this->bc[0]->set(yForm.mutable_array(), std::nullopt);
        //dolfinx::fem::petsc::set_bc<T>(yForm.mutable_array(), nullptr, this->bc);
        dolfinx::fem::petsc::set_bc<T>(this->y.vec(), this->bc, nullptr);
        yForm.scatter_fwd(); // new addition
    }

    template<typename T> LU_solver<T>::LU_solver(
        const std::shared_ptr<const varForm<T>>& aIn, 
        const std::shared_ptr<const varForm<T>>& LIn, 
        const std::vector<std::shared_ptr<const dolfinxDirichletBC<T>>>& bcIn,  
        const std::shared_ptr<dolfinxFunction<T>>& solFuncIn,
        const bool& useMumps):
        linear_solver<T>(aIn, LIn, bcIn, solFuncIn){

            dolfinx::la::petsc::options::set("ksp_type", "preonly");
            dolfinx::la::petsc::options::set("pc_type", "lu");
            if (useMumps == true){
                dolfinx::la::petsc::options::set("pc_factor_mat_solver_type", "mumps");
            }
            this->linearSolver.set_from_options();
            this->linearSolver.set_operator(this->A.mat());
    }

    template<typename T> void LU_solver<T>::solve(){
        this->linearSolver.solve(this->x.vec(), this->y.vec());
    }

    template class linear_solver<PetscScalar>;
    template class LU_solver<PetscScalar>;
}