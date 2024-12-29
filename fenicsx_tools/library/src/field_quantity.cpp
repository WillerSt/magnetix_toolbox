#include <magnetics_toolbox/field_quantity.h>

namespace mag_tools{
    template<typename T> linear_superposition<T>::linear_superposition(const std::vector<std::shared_ptr<dolfinxFunction<T>>>& inputFunctions,  std::vector<double> scaleIn):
        field_quantity<T>(inputFunctions[0]->function_space()), inFcts(inputFunctions),  
        scale(scaleIn), outVec(dolfinxVector(la::petsc::create_vector_wrap(*this->outFct->x()), false)){

    for (std::size_t i = 0; i < inFcts.size(); i++){
        dofVecs.push_back(dolfinxVector(la::petsc::create_vector_wrap(*inFcts[i]->x()), false));
    }
    update_result();

    }

    template<typename T> void linear_superposition<T>::update_result(){
        this->outFct->x()->set(0.0);
        for (std::size_t i = 0; i < this->inFcts.size(); i++){
            VecAXPY(outVec.vec(), scale[i], this->dofVecs[i].vec());
        }
        VecGhostUpdateBegin(outVec.vec(), INSERT_VALUES, SCATTER_FORWARD);
        VecGhostUpdateEnd(outVec.vec(), INSERT_VALUES, SCATTER_FORWARD);
            
    }

    template<typename T> surf_curl_evaluation<T>::surf_curl_evaluation(const std::shared_ptr<dolfinxFunction<T>>& solIn, const std::shared_ptr<const dolfinxFS<T>>& quadFSIn, const T& scaleIn):
            field_quantity<T>(std::make_shared<dolfinxFunction<T>>(quadFSIn)), quadFS(this->outFct->function_space()), sol(solIn), scale(scaleIn){
                this->initialize_variables();
                this->update_quad_point_values();
            }

    template<typename T> surf_curl_evaluation<T>::surf_curl_evaluation(const std::shared_ptr<dolfinxFunction<T>>& solIn, const std::shared_ptr<dolfinxFunction<T>>& quadFuncIn, const T& scaleIn):
             field_quantity<T>(quadFuncIn), quadFS(quadFuncIn->function_space()), sol(solIn), scale(scaleIn){
                this->initialize_variables();
                this->update_quad_point_values();
            }

    template<typename T> void surf_curl_evaluation<T>::update_result(){
        this->update_quad_point_values();
    }

    #ifdef USE_QUICK_VERSION
    template<typename T> void surf_curl_evaluation<T>::initialize_variables(){
        // assemble projection matrix
        MatZeroEntries(curlProjMat.mat());
        fem::assemble_matrix(la::petsc::Matrix::set_block_fn(curlProjMat.mat(), ADD_VALUES), *projForm, {});
        MatAssemblyBegin(curlProjMat.mat(), MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(curlProjMat.mat(), MAT_FINAL_ASSEMBLY);

        // assemble evaluation matrix and get diagonal
        std::shared_ptr<varForm<T>> evalForm = std::make_shared<varForm<T>>(dolfinx::fem::create_form<T>(
            *form_curl_evaluation_a, {quadFS, quadFS}, {{"dummy",std::shared_ptr<const dolfinxFunction<T>>()}}, {}, {}, {}));
        dolfinxMatrix curlEvalMat = dolfinx::la::petsc::Matrix(dolfinx::fem::petsc::create_matrix(*evalForm), false);
        MatZeroEntries(curlEvalMat.mat());
        fem::assemble_matrix(la::petsc::Matrix::set_block_fn(curlEvalMat.mat(), ADD_VALUES), *evalForm, {});
        MatAssemblyBegin(curlEvalMat.mat(), MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(curlEvalMat.mat(), MAT_FINAL_ASSEMBLY);                                   
        
        MatGetDiagonal(curlEvalMat.mat(), this->quadDiagVec.vec());
        for (std::size_t i = 0; i< quadDiag.x()->mutable_array().size(); i++){
            quadDiag.x()->mutable_array()[i]=1/quadDiag.x()->mutable_array()[i];
        }
    }
    #else
    template<typename T> void surf_curl_evaluation<T>::initialize_variables(){
        solver.assemble_matrix();
    }
    #endif

    #ifdef USE_QUICK_VERSION
    template<typename T> void surf_curl_evaluation<T>::update_quad_point_values(){
        MatMult(curlProjMat.mat(), solVec.vec(), curlQuadVec.vec());
        VecPointwiseMult(curlQuadVec.vec(), quadDiagVec.vec(), curlQuadVec.vec());
    }
    #else
    template<typename T> void surf_curl_evaluation<T>::update_quad_point_values(){
        solver.assemble_rhs();
        solver.solve();
    }
    #endif

    template class linear_superposition<PetscScalar>;
    template class surf_curl_evaluation<PetscScalar>;
}