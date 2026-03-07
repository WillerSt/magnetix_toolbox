// SPDX-FileCopyrightText: 2025 Stephan Willerich
// SPDX-License-Identifier: MIT License

#pragma once
#include "magnetics_toolbox/mag_tools_basic.h"
#include <dolfinx_mpc/ContactConstraint.h>
#include <dolfinx_mpc/MultiPointConstraint.h>
#include <dolfinx_mpc/assemble_matrix.h>
#include <dolfinx_mpc/utils.h>
#include <dolfinx_mpc/lifting.h>
#include <dolfinx_mpc/assemble_vector.h>

#include <dolfinx.h>
#include <dolfinx/fem/petsc.h>

namespace mag_tools{
    //template <typename T> 
    class constrained_solver{

        using T = PetscScalar;
        protected:
        const std::shared_ptr<const varForm<T>> a;
        const std::shared_ptr<const varForm<T>> L;
        const std::shared_ptr<std::vector<dolfinx::fem::DirichletBC<T>>> bc;   

        const std::shared_ptr<dolfinxFunction<T>> solFunc;

        const std::shared_ptr<const dolfinx::mesh::MeshTags<int32_t>> meshMarkers;

        

        dolfinxRHS<T> yForm = dolfinx::la::Vector<T>(L->function_spaces()[0]->dofmap()->index_map,
                    L->function_spaces()[0]->dofmap()->index_map_bs());
        
        dolfinxLinearSolver linearSolver = dolfinxLinearSolver(MPI_COMM_WORLD);   
        
        dolfinxMatrix A = dolfinx::la::petsc::Matrix(dolfinx::fem::petsc::create_matrix(*a), false);
        dolfinxVector x = dolfinxVector(dolfinx::la::petsc::create_vector_wrap(*solFunc->x()), false); 
        dolfinxVector y = dolfinxVector(dolfinx::la::petsc::create_vector_wrap(yForm), false); 

        bool mpcPresent = false;
        int idxMaster;
        int idxSlave;
        double tol;

        std::shared_ptr<dolfinx_mpc::mpc_data<T>> mpcInput = nullptr;
        std::shared_ptr<dolfinx_mpc::MultiPointConstraint<T,U<T>>> mpc = nullptr;     

        public:

        constrained_solver(
            const std::shared_ptr<const varForm<T>>& aIn, 
            const std::shared_ptr<const varForm<T>>& LIn, 
            const std::shared_ptr<std::vector<dolfinx::fem::DirichletBC<T>>>& bcIn,  
            const std::shared_ptr<dolfinxFunction<T>>& solFuncIn,
            const std::shared_ptr<const dolfinx::mesh::MeshTags<int32_t>> meshMarkersIn):
            a(aIn), L(LIn), bc(bcIn), solFunc(solFuncIn), meshMarkers(meshMarkersIn)
            {
            dolfinx::la::petsc::options::set("ksp_type", "preonly");
            dolfinx::la::petsc::options::set("pc_type", "lu");

            this->linearSolver.set_from_options();
            //this->set_mpc(6,7);
            
        }

        void set_mpc(const int& idxMasterIn, const int& idxSlaveIn, const double& tolIn){
            idxMaster = idxMasterIn;
            idxSlave = idxSlaveIn;
            tol = tolIn;
            mpcPresent = true;
            update_mpc();
        }


        void update_mpc(){
            if (this->mpcPresent){
                //this->mpc = std::make_shared<dolfinx_mpc::MultiPointConstraint<T,U<T>>> (dolfinx_mpc::MultiPointConstraint<T,U<T>>((a->function_spaces()[0]),{},{},{},{},{}));
                this->mpcInput = std::make_shared<dolfinx_mpc::mpc_data<T>>(dolfinx_mpc::create_contact_inelastic_condition<T,U<T>>(*solFunc->function_space(), *meshMarkers, idxSlave, idxMaster, tol));
                this->mpc = std::make_shared<dolfinx_mpc::MultiPointConstraint<T,U<T>>> ( dolfinx_mpc::MultiPointConstraint<T,U<T>>(
                    solFunc->function_space(),
                    mpcInput->slaves,
                    mpcInput->masters,
                    mpcInput->coeffs,
                    mpcInput->owners,
                    mpcInput->offsets));
            }
        }



        void assemble_matrix(){
            if (this->mpcPresent){
                A = dolfinx_mpc::create_matrix(*a, mpc);
            }         

            // build vector of reference_wrappers for the DirichletBCs to match assembler overload



            
            
            MatZeroEntries(A.mat());  

            if (this->mpcPresent){    
                using BC_element_t = dolfinx::fem::DirichletBC<T>;
                std::vector<std::shared_ptr<const BC_element_t>> bc_refs;
                bc_refs.reserve(bc->size());
                for (auto &d : *bc)
                    bc_refs.emplace_back(std::make_shared<const BC_element_t>(d));  
                dolfinx_mpc::assemble_matrix(
                    dolfinx::la::petsc::Matrix::set_block_fn(A.mat(), ADD_VALUES),
                    dolfinx::la::petsc::Matrix::set_fn(A.mat(), ADD_VALUES),*a, 
                    mpc, 
                    mpc,
                    bc_refs);
            }
            else{
                using BC_element_t = std::remove_reference_t<decltype((*bc)[0])>;
                std::vector<std::reference_wrapper<const BC_element_t>> bc_refs;
                bc_refs.reserve(bc->size());
                for (auto &d : *bc)
                    bc_refs.emplace_back(std::cref(d));
                dolfinx::fem::assemble_matrix(dolfinx::la::petsc::Matrix::set_block_fn(A.mat(), ADD_VALUES),
                        *a, bc_refs);
            }

            MatAssemblyBegin(A.mat(), MAT_FLUSH_ASSEMBLY);
            MatAssemblyEnd(A.mat(), MAT_FLUSH_ASSEMBLY);

            using BC_element_t = std::remove_reference_t<decltype((*bc)[0])>;
                std::vector<std::reference_wrapper<const BC_element_t>> bc_refs;
                bc_refs.reserve(bc->size());
                for (auto &d : *bc)
                    bc_refs.emplace_back(std::cref(d));
                    
            dolfinx::fem::set_diagonal<T>(dolfinx::la::petsc::Matrix::set_fn(A.mat(), INSERT_VALUES), 
                            *(a->function_spaces()[0]),
                            bc_refs);
            MatAssemblyBegin(A.mat(), MAT_FINAL_ASSEMBLY);
            MatAssemblyEnd(A.mat(), MAT_FINAL_ASSEMBLY);
        }
        void assemble_rhs(){
            this->assemble_rhs_internally(this->L);
        }

        double calc_orig_rhs_norm(const std::shared_ptr<const varForm<T>>& zeroL){
            // assemble_rhs() has to be called externally to get back to old state
            assemble_rhs_internally(zeroL);
            return get_rhs_norm();         
        }


        double get_rhs_norm() const{
            return dolfinx::la::norm(yForm);
        }

        void solve(){
            this->linearSolver.set_operator(this->A.mat());
            this->linearSolver.solve(this->x.vec(), this->y.vec());

            if (mpcPresent){
                mpc->backsubstitution(this->solFunc->x()->array());
            }

        }
        //inline virtual ~linear_solver(){};

        protected:
        void assemble_rhs_internally(const std::shared_ptr<const varForm<T>>& formIn){
            std::ranges::fill(yForm.array(), 0.0);

            // build vector of reference_wrappers for the DirichletBCs to match assembler overload

            if (this->mpcPresent == true){
                using BC_element_t = dolfinx::fem::DirichletBC<T>;
                std::vector<std::shared_ptr<const BC_element_t>> bc_refs;
                bc_refs.reserve(bc->size());
                for (auto &d : *bc)
                    bc_refs.emplace_back(std::make_shared<const BC_element_t>(d));  
                dolfinx_mpc::assemble_vector(yForm.array(), *formIn, mpc);  
                dolfinx_mpc::apply_lifting(yForm.array(), {this->a}, {bc_refs}, {}, T(1.0), mpc);
            }
            else{

                using BC_element_t = std::remove_reference_t<decltype((*bc)[0])>;
                std::vector<std::reference_wrapper<const BC_element_t>> bc_refs = {};
                bc_refs.reserve(bc->size());
                for (auto &d : *bc){
                    bc_refs.emplace_back(std::cref(d));
                }
                
                std::ranges::fill(yForm.array(), 0.0);
                dolfinx::fem::assemble_vector(yForm.array(), *formIn);

                dolfinx::fem::apply_lifting(yForm.array(), {*this->a}, {bc_refs}, std::vector<std::span<const T>>{}, T(1.0));
                yForm.scatter_rev(std::plus<T>());

            }          
            
            for (auto &d : *bc)
                d.set(yForm.array(), std::nullopt, T(1.0));
                
            yForm.scatter_fwd(); // new addition   
        }


    };

}