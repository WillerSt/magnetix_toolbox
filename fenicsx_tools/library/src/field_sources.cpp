// SPDX-FileCopyrightText: 2025 Stephan Willerich
// SPDX-License-Identifier: MIT License

#include <magnetics_toolbox/field_sources.h>
#include <magnetics_toolbox/constants.h>

namespace mag_tools::src{
    sine_variation::sine_variation(const double& ampIn, const double& fIn, const double& phaseIn, const double& offsetIn)
        :amp(ampIn), f(fIn), phase(phaseIn), offset(offsetIn){}

    double sine_variation::calc_current(const double& t) const{
        return this->amp * std::sin(2*constants::pi * this->f * t + this->phase) + this->offset;
    }


    template<typename T> const_scalar_source<T>::const_scalar_source(const std::shared_ptr<dolfinxFunction<T>>& quadFuncIn, const std::shared_ptr<const meshTags>& meshMarkers , const std::vector<double>& scaleIn,
            const std::vector<int>& tagIdxIn, const std::shared_ptr<const time_variation>& curFuncIn):quadFunc(quadFuncIn), scaleFactors(scaleIn), tagIdx(tagIdxIn), curFunc(curFuncIn){
                
                for (std::size_t i = 0; i < this->tagIdx.size(); i++){
                    coupledDofs.push_back(quadrature_dof_coupler<T>(quadFunc, meshMarkers, this->tagIdx[i], 1, 1));
                }
    }

     template<typename T> void const_scalar_source<T>::update_source(const double& t){
        this->curVal = this->curFunc->calc_current(t);
        this->set_dofs_to_curVal();
    }

     template<typename T> double const_scalar_source<T>::get_excitation_value() const{
        return this->curVal;
    }

    template<typename T> void const_scalar_source<T>::set_excitation_unity(){
        this->curVal = 1.0;
        this->set_dofs_to_curVal();
    }

    template<typename T> void const_scalar_source<T>::set_excitation_zero(){
        this->curVal = 0.0;
        this->set_dofs_to_curVal();
    }


    template<typename T> void const_scalar_source<T>::set_dofs_to_curVal(){
            for (std::size_t i = 0; i < this->coupledDofs.size(); i++){
            this->coupledDofs[i].set_all_entries(this->scaleFactors[i] * this->curVal);
        }
    }

    template class const_scalar_source<PetscScalar>;
}