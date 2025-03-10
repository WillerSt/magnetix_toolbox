// SPDX-FileCopyrightText: 2025 Stephan Willerich
// SPDX-License-Identifier: MIT License

#include <magnetics_toolbox/function_processing.h>

#include <quadVec_dgVec_rotated.h>

namespace mag_tools{

    template <typename T> std::shared_ptr<const dolfinxFunction<T>> create_const_function(const std::shared_ptr<const dolfinxFS<T>>& funcSpace, const T& value){
        auto cFunc  = std::make_shared<dolfinxFunction<T>>(funcSpace);
        cFunc->x()->set(value);
        return cFunc;
    }

    template <typename T> dg_projection<T>::dg_projection(const std::shared_ptr<dolfinxMesh<T>>& meshIn, const int& dimIn):
        dgVecFS(initialze_function_space(meshIn, dimIn)),
        a(initialize_form(dimIn)), dim(dimIn){

        MatZeroEntries(A.mat());
        dolfinx::fem::assemble_matrix(la::petsc::Matrix::set_block_fn(A.mat(), ADD_VALUES), *a, {});

        MatAssemblyBegin(A.mat(), MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(A.mat(), MAT_FINAL_ASSEMBLY);


        la::petsc::options::set("ksp_type", "preonly");
        la::petsc::options::set("pc_type", "lu");
        lu.set_from_options();

        lu.set_operator(A.mat());

    }

    template <typename T> std::shared_ptr<dolfinxFunction<T>> dg_projection<T>::project_quad_to_dg(const std::shared_ptr<const dolfinxFunction<T>>& funcIn, const std::shared_ptr<dolfinxFunction<T>>& B){
        
        std::shared_ptr<const fem::Form<T>> L;

        if (this->dim > 1){
            L = std::make_shared<const fem::Form<T>>(dolfinx::fem::create_form<T>(
                *form_quadVec_dgVec_L, {dgVecFS}, {{"sol", funcIn}}, {},{}, {}));}
        else if (this->dim ==1){
            L = std::make_shared<const fem::Form<T>>(dolfinx::fem::create_form<T>(
                *form_quadSca_dgSca_L, {dgVecFS}, {{"sol", funcIn}}, {},{}, {}));      
        
        }
        

        la::Vector<T> b(L->function_spaces()[0]->dofmap()->index_map,
                        L->function_spaces()[0]->dofmap()->index_map_bs());          

        b.set(0.0);
        fem::assemble_vector(b.mutable_array(), *L);
        b.scatter_fwd(); // seam visible -> ghost mode none gets rid of it
        la::petsc::Vector _u(la::petsc::create_vector_wrap(*B->x()), false);
        la::petsc::Vector _b(la::petsc::create_vector_wrap(b), false);
        lu.solve(_u.vec(), _b.vec());
        B->x()->scatter_fwd();
        return B;
    }

    template <typename T> std::shared_ptr<dolfinxFunction<T>> dg_projection<T>::project_quad_to_dg(const std::shared_ptr<const dolfinxFunction<T>>& funcIn, const std::shared_ptr<dolfinxFunction<T>>& B, const std::shared_ptr<dolfinxFunction<T>>& rotFunc){
        
        std::shared_ptr<const fem::Form<T>> L;

        if (this->dim > 1){
            L = std::make_shared<const fem::Form<T>>(dolfinx::fem::create_form<T>(
                *form_quadVec_dgVec_rotated_L, {dgVecFS}, {{"sol", funcIn},{"rotFunc", rotFunc}}, {},{}, {}));}
        else if (this->dim ==1){
            std::cout << "ERROR: Rotation Matrix cannot be used in conjunction with scalar quantity\n";        
        }
        

        la::Vector<T> b(L->function_spaces()[0]->dofmap()->index_map,
                        L->function_spaces()[0]->dofmap()->index_map_bs());          

        b.set(0.0);
        fem::assemble_vector(b.mutable_array(), *L);

        la::petsc::Vector _u(la::petsc::create_vector_wrap(*B->x()), false);
        la::petsc::Vector _b(la::petsc::create_vector_wrap(b), false);
        lu.solve(_u.vec(), _b.vec());
        return B;
    }


    template <typename T> std::shared_ptr<const varForm<T>> dg_projection<T>::initialize_form(const int& dim){
        if (dim == 2){
            return std::make_shared<const varForm<T>>(dolfinx::fem::create_form<T>(*form_quadVec_dgVec_a, {dgVecFS, dgVecFS}, {{"dummy",std::shared_ptr<const dolfinxFunction<T>>()}}, {}, {}, {}));
        }
        else if (dim == 1){
            return std::make_shared<const varForm<T>>(dolfinx::fem::create_form<T>(*form_quadSca_dgSca_a, {dgVecFS, dgVecFS}, {{"dummy",std::shared_ptr<const dolfinxFunction<T>>()}}, {}, {}, {}));
        }
        else{
            std::cout << "DG-Projection for dim " << dim << " is not supported";
            exit(1); 
        }
    }

    template <typename T> std::shared_ptr<dolfinxFS<T>> dg_projection<T>::initialze_function_space(const std::shared_ptr<dolfinxMesh<T>>& meshIn, const int& dim){
                // initialize function spaces and functions
        auto dg1Elemem = basix::create_element<U<T>>(
            basix::element::family::P, basix::cell::type::triangle, 1,
            basix::element::lagrange_variant::unset,
            basix::element::dpc_variant::unset, true);        
        
        if (dim == 2){
            return std::make_shared<dolfinxFS<T>>(fem::create_functionspace(meshIn, dg1Elemem, {2}));
        }
        else if (dim == 1){
            return std::make_shared<dolfinxFS<T>>(fem::create_functionspace(meshIn, dg1Elemem, {}));
        }
        else{
            std::cout << "DG-Projection for dim " << dim << " is not supported";
            exit(1); 
        }
    }

    template <typename T> std::shared_ptr<dolfinxFunction<T>> project_quad_to_dg(const std::shared_ptr<const dolfinxFunction<T>>& funcIn, const std::shared_ptr<dolfinxMesh<T>>& meshIn, const std::shared_ptr<dolfinxFunction<T>>& B){

        auto dg1Elemem = basix::create_element<U<T>>(
            basix::element::family::P, basix::cell::type::triangle, 1,
            basix::element::lagrange_variant::unset,
            basix::element::dpc_variant::unset, true);  
        auto dgVecFS = std::make_shared<dolfinxFS<T>>(fem::create_functionspace(meshIn, dg1Elemem, {2,1}));
        
        // auto B = std::make_shared<dolfinxFunction>(dgVecFS);

        auto a = std::make_shared<const varForm<T>>(dolfinx::fem::create_form<T>(
            *form_quadVec_dgVec_a, {dgVecFS, dgVecFS}, {{"dummy",std::shared_ptr<const dolfinxFunction<T>>()}}, {}, {}, {}));
        auto L = std::make_shared<const fem::Form<T>>(dolfinx::fem::create_form<T>(
            *form_quadVec_dgVec_L, {dgVecFS}, {{"sol", funcIn}}, {},{}, {}));


    
        
        auto A = dolfinx::la::petsc::Matrix(dolfinx::fem::petsc::create_matrix(*a), false);
        la::Vector<T> b(L->function_spaces()[0]->dofmap()->index_map,
                        L->function_spaces()[0]->dofmap()->index_map_bs());

        MatZeroEntries(A.mat());
        dolfinx::fem::assemble_matrix(la::petsc::Matrix::set_block_fn(A.mat(), ADD_VALUES),
                            *a, {});

        MatAssemblyBegin(A.mat(), MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(A.mat(), MAT_FINAL_ASSEMBLY);
        


        b.set(0.0);
        fem::assemble_vector(b.mutable_array(), *L);

        la::petsc::KrylovSolver lu(MPI_COMM_WORLD);
        la::petsc::options::set("ksp_type", "preonly");
        la::petsc::options::set("pc_type", "lu");
        lu.set_from_options();

        lu.set_operator(A.mat());
        la::petsc::Vector _u(la::petsc::create_vector_wrap(*B->x()), false);
        la::petsc::Vector _b(la::petsc::create_vector_wrap(b), false);
        lu.solve(_u.vec(), _b.vec());
        B->x()->scatter_fwd();
        /*
        auto Ainv = dolfinx::la::petsc::Matrix(dolfinx::fem::petsc::create_matrix(*a), false);
        MatInvertBlockDiagonalMat(A.mat(), Ainv.mat());
        MatAssemblyBegin(Ainv.mat(), MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(Ainv.mat(), MAT_FINAL_ASSEMBLY);

        auto B_alt = std::make_shared<dolfinxFunction>(dgVecFS);
        la::petsc::Vector _x(la::petsc::create_vector_wrap(*B_alt->x()), false);
        
        MatMult(Ainv.mat(), _b.vec(), _x.vec());

        for (std::size_t i = 0; i<_x.local_size();i++)
        {
            std::cout << B_alt->x()->array()[i] /B->x()->array()[i] << std::endl;
        }
        */
        return B;
    }

    template class dg_projection<PetscScalar>;
    template std::shared_ptr<dolfinxFunction<PetscScalar>> project_quad_to_dg<PetscScalar>(const std::shared_ptr<const dolfinxFunction<PetscScalar>>& funcIn, const std::shared_ptr<dolfinxMesh<PetscScalar>>& meshIn, const std::shared_ptr<dolfinxFunction<PetscScalar>>& B);
    template std::shared_ptr<const dolfinxFunction<PetscScalar>> create_const_function<PetscScalar>(const std::shared_ptr<const dolfinxFS<PetscScalar>>& funcSpace, const PetscScalar& value);
}   