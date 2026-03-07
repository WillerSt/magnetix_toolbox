// SPDX-FileCopyrightText: 2025 Stephan Willerich
// SPDX-License-Identifier: MIT License

#include <magnetics_toolbox/quadrature_dof_coupler.h>
#include <math.h>

namespace mag_tools{

    template<typename T> quadrature_dof_coupler<T>::quadrature_dof_coupler(    
            const std::shared_ptr<dolfinx::fem::Function<T>>& funcIn, 
            const std::shared_ptr<const dolfinx::mesh::MeshTags<std::int32_t>>& meshMarkersIn,
            const int& tagIn,
            const int& dofDimRowsIn,
            const int& dofDimColsIn):
            quadFunc(funcIn), meshMarkers(meshMarkersIn), tag(tagIn), dofDimRows(dofDimRowsIn), dofDimCols(dofDimColsIn),
                refVector(std::make_shared<refVecType>()), 
                idxVector(std::make_shared<idxVecType>()){
            
            std::vector<int32_t>  cellIndices = meshMarkers->find(tag);
            this->initialize_coupling(cellIndices);   

    }


    template<typename T> quadrature_dof_coupler<T>::quadrature_dof_coupler(    
            const std::shared_ptr<dolfinx::fem::Function<T>>& funcIn,
            const size_t& dofDimRowsIn,
            const size_t& dofDimColsIn):
            quadFunc(funcIn), meshMarkers(nullptr), tag(-1), dofDimRows(dofDimRowsIn), dofDimCols(dofDimColsIn),
                refVector(std::make_shared<refVecType>()), 
                idxVector(std::make_shared<idxVecType>()){

            std::vector<int32_t>  cellIndices;         
            
            for (int i = 0; i < this->quadFunc->function_space()->mesh()->topology()->index_map(this->quadFunc->function_space()->mesh()->topology()->dim())->size_local(); i++){
                cellIndices.push_back(i);
            }

            this->initialize_coupling(cellIndices); 

    }

    template<typename T> int quadrature_dof_coupler<T>::report_number_of_dofs(bool display) const{
        int noOfDofs; 
        if (this->refVector->size() >0){
            noOfDofs =  this->refVector->size() * (*this->refVector)[0].size() * (*this->refVector)[0][0].size();
        }
        else{
            noOfDofs = 0;
        }
        if (display==true){
                std::cout << "No. of dofs: "<< noOfDofs <<"\n";
        }
        return noOfDofs;
    }


    
    template<typename T>  void quadrature_dof_coupler<T>::perform_checks()const{

    }

    template<typename T>  void quadrature_dof_coupler<T>::coordinate_check() const{
        for (unsigned int dof = 0; dof < this->idxVector->size(); dof++){
            for (size_t j = 0; j < this->dofDimCols; j++){
                for (size_t i = 0; i < this->dofDimRows; i++){
                    
                }

            }
        }
    }

    template<typename T> void quadrature_dof_coupler<T>::set_all_entries_zero() const{
        this->set_all_entries(0.0);
    }

    template<typename T> void quadrature_dof_coupler<T>::set_all_entries(const T& valIn) const{
        for (auto& dof: *(this->refVector)){
            for (size_t i = 0; i < this->dofDimRows; i++){
                for(size_t j = 0; j<this->dofDimCols; j++){
                    *(dof[i][j]) = valIn;
                }
            }
        }
    }

    template<typename T> void quadrature_dof_coupler<T>::set_diagonal(const T& valIn) const{
        if (this->dofDimCols == this->dofDimRows)
        {
            for (auto& dof: *(this->refVector)){
                for (size_t i = 0; i < this->dofDimRows; i++){
                    for(size_t j = 0; j<this->dofDimCols; j++){
                        if (i==j){
                            *(dof[i][j]) = valIn;
                        }
                    }
                }
            }
        }
    }
    
    template<typename T> void quadrature_dof_coupler<T>::set_all_indexed_values(const size_t& row, const size_t& col, const T& valIn) const{
        if ((col < this->dofDimCols) &&  (row < this->dofDimRows))
        {
            for (auto& dof: *(this->refVector)){
                *(dof[row][col]) = valIn;
            }
        }
        else{
            std::cout << "unable to set values, row: " << row << " dim " << this->dofDimRows << ", col " << col << " dim " << this->dofDimCols <<std::endl;
        }

    }

    template<typename T> std::shared_ptr<const typename quadrature_dof_coupler<T>::refVecType> quadrature_dof_coupler<T>::get_ref_vec() const{
        return this->refVector;
    }

    template<typename T> std::shared_ptr<const typename quadrature_dof_coupler<T>::idxVecType> quadrature_dof_coupler<T>::get_idx_vec() const{
        return this->idxVector;

    }

    template<typename T> void quadrature_dof_coupler<T>::initialize_coupling(const std::vector<int32_t>&  cellIndices){
        /*
            1. loop over all cells
                2. loop over components
                    components numbered according to scheme:
                    -> 1 Dim x
                    -> 2 Dim    Vector:
                                        c0
                                        c1
                                Tensor / Matrix:
                                        c0  c2
                                        c1  c3
                    -> 3 Dim    Vector:
                                        c0
                                        c1
                                        c2
                                Tensor / Matrix:
                                        c0  c3  c6
                                        c1  c4  c7
                                        c2  c5  c8
                -> result: indices of compoenents in underlying PETSc vector in  variable componentDofs
                        Point 0     Point 1     ....    Point n
                        c0          c0          ....    c0
                        c1          c1          ....    c1
                        ...
                        cn          cn          ....    cn    
                
                
                3. loop over points within cell
                    reorder variable componentDofs created in step 2
                    -> pointDofs: Pointer to dofs associated with a point order as vector / matrix
                    -> idxMat: indices dofs associated with a point order as vector / matrix
            
             create data structures containing pointers and indices of all dofs ordered as vector / matrix
            -> refVector pointer to values in PETSc vector
            -> idx Vector indices of values in PETSc vector
                -> access according to refVector[pointIdx][rowIdx][colIdx]
                    
                    

        */

        auto shape = this->quadFunc->function_space()->element()->value_shape();
        
        bool matCheck = shape.size()>1;
        auto bs = this->quadFunc->function_space()->dofmap()->bs();
        
        size_t noC;
        if (shape.size()>0){
            noC = shape[0];
        }
        else{
            noC = 1;
        }
        
        if (matCheck){
            noC = noC * shape[1];
        }
        
        
        for (auto& cell:cellIndices){
            auto cellDofs = this->quadFunc->function_space()->dofmap()->cell_dofs(cell);

            if (cellDofs.size()*bs%noC != 0){
               std::cout << "No of cell dofs " << cellDofs.size() << " not a multiple of no of expected components " << noC << std::endl;
            }


            size_t noU = cellDofs.size();

            for (size_t uIdx = 0; uIdx<noU; uIdx++){
                std::vector<std::vector<int>> idxMat = {};
                std::vector<std::vector<T*>> dofMat ={};
            
                for (size_t i = 0; i<dofDimRows; i++){
                        std::vector<int> rowIdxTemp = {};
                        std::vector<T*> rowPtrTemp  = {};
                        
                    for (size_t j = 0; j<dofDimCols; j++){

                            rowIdxTemp.push_back(-1);
                            rowPtrTemp.push_back(nullptr);
                        }
                    idxMat.push_back(rowIdxTemp);
                    dofMat.push_back(rowPtrTemp);
                }
                for (size_t i = 0; i<dofDimRows; i++){
                    for (size_t j = 0; j<dofDimCols; j++){
                        idxMat[i][j] = cellDofs[uIdx]*bs+dofDimCols*i+j;
                    }
                }
                
                
                for (size_t i = 0; i<dofDimRows; i++){
                    for (size_t j = 0; j<dofDimCols; j++){
                        dofMat[i][j] = &(this->quadFunc->x()->array()[idxMat[i][j]]);
                    }
                }

                this->idxVector->push_back(idxMat);
                this->refVector->push_back(dofMat);
                
            }
            
        }
    }


    template class quadrature_dof_coupler<PetscScalar>;
}