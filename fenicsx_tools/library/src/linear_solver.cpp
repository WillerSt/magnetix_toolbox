// SPDX-FileCopyrightText: 2025 Stephan Willerich
// SPDX-License-Identifier: MIT License

#include <magnetics_toolbox/linear_solver.h>
#include <dolfinx/fem/petsc.h>

namespace mag_tools{
    template<typename T> void linear_solver<T>::assemble_matrix(){
        MatZeroEntries(A.mat());

        // build vector of reference_wrappers for the DirichletBCs to match assembler overload
        using BC_element_t = std::remove_reference_t<decltype((*bc)[0])>;
        std::vector<std::reference_wrapper<const BC_element_t>> bc_refs;
        bc_refs.reserve(bc->size());
        for (auto &d : *bc)
            bc_refs.emplace_back(std::cref(d));

        dolfinx::fem::assemble_matrix(dolfinx::la::petsc::Matrix::set_block_fn(A.mat(), ADD_VALUES),
                        *a, bc_refs);
        MatAssemblyBegin(A.mat(), MAT_FLUSH_ASSEMBLY);
        MatAssemblyEnd(A.mat(), MAT_FLUSH_ASSEMBLY);
            
        dolfinx::fem::set_diagonal<T>(dolfinx::la::petsc::Matrix::set_fn(A.mat(), INSERT_VALUES), 
                        *(a->function_spaces()[0]),
                        bc_refs);
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
        std::ranges::fill(yForm.array(), 0.0);
        dolfinx::fem::assemble_vector(yForm.array(), *formIn);

        // build vector of reference_wrappers for the DirichletBCs to match assembler overload
        using BC_element_t = std::remove_reference_t<decltype((*bc)[0])>;
        std::vector<std::reference_wrapper<const BC_element_t>> bc_refs = {};
        bc_refs.reserve(bc->size());
        for (auto &d : *bc){
            bc_refs.emplace_back(std::cref(d));
        }

        dolfinx::fem::apply_lifting(yForm.array(), {*this->a}, {bc_refs}, std::vector<std::span<const T>>{}, T(1.0));
        yForm.scatter_rev(std::plus<T>());

        for (auto &d : *bc)
            d.set(yForm.array(), std::nullopt, T(1.0));
            
        yForm.scatter_fwd(); // new addition
    }

    template<typename T> LU_solver<T>::LU_solver(
        const std::shared_ptr<const varForm<T>>& aIn, 
        const std::shared_ptr<const varForm<T>>& LIn, 
        const std::shared_ptr<std::vector<dolfinx::fem::DirichletBC<T>>>& bcIn,  
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