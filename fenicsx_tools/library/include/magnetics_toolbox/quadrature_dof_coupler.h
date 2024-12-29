#pragma once
#include <magnetics_toolbox/mag_tools_basic.h>


namespace mag_tools{
template<typename T> class quadrature_dof_coupler{

public:

using refVecType = std::vector<std::vector<std::vector<T*>>>;
using idxVecType = std::vector<std::vector<std::vector<int>>>;

public:



const std::shared_ptr<dolfinxFunction<T>> quadFunc;
const std::shared_ptr<const mesh::MeshTags<std::int32_t>> meshMarkers;

const int tag;
const size_t dofDimRows;
const size_t dofDimCols;

//const int cSize;

const std::shared_ptr<refVecType> refVector;
const std::shared_ptr<idxVecType> idxVector;

public:
    quadrature_dof_coupler(    
        const std::shared_ptr<dolfinx::fem::Function<T>>& funcIn, 
        const std::shared_ptr<const mesh::MeshTags<std::int32_t>>& meshMarkersIn,
        const int& tagIn,
        const int& dofDimRowsIn,
        const int& dofDimColsIn);

    quadrature_dof_coupler(    
        const std::shared_ptr<dolfinx::fem::Function<T>>& funcIn,
        const size_t& dofDimRowsIn,
        const size_t& dofDimColsIn);

    int report_number_of_dofs(bool display = false) const;

    void perform_checks()const;

    void coordinate_check() const;

    void set_all_entries_zero() const;

    void set_all_entries(const T& valIn) const;

    void set_diagonal(const T& valIn) const;
    
    void set_all_indexed_values(const size_t& row, const size_t& col, const T& valIn) const;

    std::shared_ptr<const refVecType> get_ref_vec() const;
    std::shared_ptr<const idxVecType> get_idx_vec() const;
private:
    void initialize_coupling(const std::vector<int32_t>&  cellIndices);

    //int determine_cSize();

};
}