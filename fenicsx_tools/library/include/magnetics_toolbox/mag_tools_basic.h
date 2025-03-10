// SPDX-FileCopyrightText: 2025 Stephan Willerich
// SPDX-License-Identifier: MIT License


#pragma once
#include <dolfinx.h>
#include <basix/quadrature.h>
#include <basix/e-lagrange.h>
#include <Eigen/Core>

namespace mag_tools{
    template <typename T> using dolfinxFunction = dolfinx::fem::Function<T>;
    template <typename T> using U = typename dolfinx::scalar_value_type_t<T>;
    template <typename T> using varForm = dolfinx::fem::Form<T,U<T>>;
    template <typename T> using dolfinxFS = dolfinx::fem::FunctionSpace<U<T>>;
    template <typename T> using dolfinxFunction = dolfinx::fem::Function<T>;
    template <typename T> using dolfinxMesh = dolfinx::mesh::Mesh<U<T>>;
    template <typename T> using dolfinxRHS = dolfinx::la::Vector<T>;
    template <typename T> using dolfinxConstant = dolfinx::fem::Constant<T>;
    template <typename T> using dolfinxDirichletBC = dolfinx::fem::DirichletBC<T>;

    using dolfinxVector = dolfinx::la::petsc::Vector;
    using dolfinxLinearSolver = dolfinx::la::petsc::KrylovSolver;
	using dolfinxMatrix = dolfinx::la::petsc::Matrix;


}

