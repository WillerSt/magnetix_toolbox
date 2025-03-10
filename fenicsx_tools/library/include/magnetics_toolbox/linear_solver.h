// SPDX-FileCopyrightText: 2025 Stephan Willerich
// SPDX-License-Identifier: MIT License

#pragma once
#include <magnetics_toolbox/mag_tools_basic.h>
#include <dolfinx/fem/petsc.h>

namespace mag_tools{
    template <typename T> class linear_solver{

        protected:
        const std::shared_ptr<const varForm<T>> a;
        const std::shared_ptr<const varForm<T>> L;
        const std::vector<std::shared_ptr<const dolfinxDirichletBC<T>>> bc;

        const std::shared_ptr<dolfinxFunction<T>> solFunc;

        dolfinxRHS<T> yForm = dolfinx::la::Vector<T>(L->function_spaces()[0]->dofmap()->index_map,
                    L->function_spaces()[0]->dofmap()->index_map_bs());

        dolfinxMatrix A = dolfinx::la::petsc::Matrix(dolfinx::fem::petsc::create_matrix(*a), false);
        dolfinxVector x = dolfinxVector(la::petsc::create_vector_wrap(*solFunc->x()), false); 
        dolfinxVector y = dolfinxVector(la::petsc::create_vector_wrap(yForm), false); 

        dolfinxLinearSolver linearSolver = dolfinxLinearSolver(MPI_COMM_WORLD);

        public:

        linear_solver(
            const std::shared_ptr<const varForm<T>>& aIn, 
            const std::shared_ptr<const varForm<T>>& LIn, 
            const std::vector<std::shared_ptr<const dolfinxDirichletBC<T>>>& bcIn,  
            const std::shared_ptr<dolfinxFunction<T>>& solFuncIn):
            a(aIn), L(LIn), bc(bcIn), solFunc(solFuncIn)
        {}

        void assemble_matrix();
        void assemble_rhs();

        double calc_orig_rhs_norm(const std::shared_ptr<const varForm<T>>& zeroL);


        double get_rhs_norm() const;

        virtual void solve() = 0;
        inline virtual ~linear_solver(){};

        protected:
        void assemble_rhs_internally(const std::shared_ptr<const varForm<T>>& formIn);


    };
    template <typename T> class LU_solver  : public linear_solver<T>{

        public:

        LU_solver(
            const std::shared_ptr<const varForm<T>>& aIn, 
            const std::shared_ptr<const varForm<T>>& LIn, 
            const std::vector<std::shared_ptr<const dolfinxDirichletBC<T>>>& bcIn,  
            const std::shared_ptr<dolfinxFunction<T>>& solFuncIn,
            const bool& useMumps = false);

        void solve();    
    };

}