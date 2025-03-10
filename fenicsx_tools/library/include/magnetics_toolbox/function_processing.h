// SPDX-FileCopyrightText: 2025 Stephan Willerich
// SPDX-License-Identifier: MIT License

#pragma once

#include <magnetics_toolbox/mag_tools_basic.h>
#include <curlcg2_dgVec.h>
#include <magnetics_toolbox/linear_solver.h>
#include <quadVec_dgVec.h>
#include <quadSca_dgSca.h>

namespace mag_tools{


    template <typename T> std::shared_ptr<const dolfinxFunction<T>> create_const_function(const std::shared_ptr<const dolfinxFS<T>>& funcSpace, const T& value);

    template <typename T> class dg_projection{
        private:
        const std::shared_ptr<dolfinxFS<T>> dgVecFS;
        const std::shared_ptr<const varForm<T>> a;
        dolfinx::la::petsc::Matrix A = dolfinx::la::petsc::Matrix(dolfinx::fem::petsc::create_matrix(*a), false);
        la::petsc::KrylovSolver lu = la::petsc::KrylovSolver(MPI_COMM_WORLD);
        const int dim;

        public:
        dg_projection(const std::shared_ptr<dolfinxMesh<T>>& meshIn, const int& dim);

        public:
        std::shared_ptr<dolfinxFunction<T>> project_quad_to_dg(const std::shared_ptr<const dolfinxFunction<T>>& funcIn, const std::shared_ptr<dolfinxFunction<T>>& B);
        std::shared_ptr<dolfinxFunction<T>> project_quad_to_dg(const std::shared_ptr<const dolfinxFunction<T>>& funcIn, const std::shared_ptr<dolfinxFunction<T>>& B, const std::shared_ptr<dolfinxFunction<T>>& rotFunc);
        private:
         std::shared_ptr<const varForm<T>> initialize_form(const int& dim);
         std::shared_ptr<dolfinxFS<T>> initialze_function_space(const std::shared_ptr<dolfinxMesh<T>>& meshIn, const int& dim);
    };

    template <typename T> std::shared_ptr<dolfinxFunction<T>> project_quad_to_dg(const std::shared_ptr<const dolfinxFunction<T>>& funcIn, const std::shared_ptr<dolfinxMesh<T>>& meshIn, const std::shared_ptr<dolfinxFunction<T>>& B);
      
}