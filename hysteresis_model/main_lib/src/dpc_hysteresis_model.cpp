// SPDX-FileCopyrightText: 2025 Stephan Willerich
// SPDX-License-Identifier: MIT License

#include <Hysteresis_Operator/dpc_hysteresis_model.h>
#include <iostream>


template <int dimension, typename CS>
DPC_hysteresis_model<dimension, CS>::DPC_hysteresis_model(const std::shared_ptr<const malloc_grid<dimension, CS>>& malloc_gridIn)
: discreteGrid(malloc_gridIn), BSat(malloc_gridIn->Bsat), hysterons(this->initialize_hysterons()),
  saturationSurfaces(discreteGrid->get_pointers_to_saturation_surfaces())
{

	for (std::size_t i = 0; i < hysterons.size(); i++ )
	{
		if (i != hysterons[i].get_associated_surfaceID())
		{
			std::cout << "ERROR: Index of hysteron " << i << " does not match ID " <<  hysterons[i].get_associated_surfaceID() << std::endl;
		}
	}
}

template DPC_hysteresis_model<2, CS_sphere<2>>::DPC_hysteresis_model(const std::shared_ptr<const malloc_grid<2, CS_sphere<2>>>& malloc_gridIn);

template <int dimension, typename CS>
std::vector<elementary_hysteron<dimension, CS>> DPC_hysteresis_model<dimension, CS>::initialize_hysterons()
{
	std::vector<elementary_hysteron<dimension, CS>> result;
	for (unsigned int i = 0; i < (*discreteGrid).numberOfSurfaces; i++)
	{
			result.push_back(elementary_hysteron<dimension, CS>(
					this->Hnew,
					this->discreteGrid->get_pointer_to_surface(i)));
			this->numberOfHysterons += 1;
			if (i != result[i].get_associated_surfaceID())
			{
				std::cout << "WARNING: Indices of hysterons are not equal to IDs of critical surfaces \n";
			}
	}

	return result;
}


template <int dimension, typename CS>
void DPC_hysteresis_model<dimension, CS>::calculate_evolution(const EigenType& HnewIn)
{
	Hnew = HnewIn;
	Hdelta = Hnew -Hold;
	value = EigenType::Zero();

	for (int i = 0; i < DPC_hysteresis_model::numberOfHysterons; i++)
	{
		value += hysterons[i].update_value(Hnew,  Hdelta, Hold);
	}

	this->value *= this->BSat;
	value +=  this->discreteGrid->muRmu0 * HnewIn;
}

template<int dimension, typename CS>
void DPC_hysteresis_model<dimension, CS>::calculate_evolution(const EigenType& HnewIn, const std::vector<int>& FrozenIDs, const EigenType& Jconst)
{
	this->Hnew = HnewIn;
	this->Hdelta = this->Hnew -this->Hold;
	this->value = EigenType::Zero();
	int calculationCounter = 0;

	for (unsigned int i = 0; i < FrozenIDs.size(); i++)
	{
		if (!hysterons[FrozenIDs[i]].frozen)
		{
			calculationCounter +=1;
			value += hysterons[FrozenIDs[i]].update_frozen_value(Hnew,  Hdelta, Hold);
		}
		else
		{
			value += hysterons[FrozenIDs[i]].oldValue;
		}
	}

	value *= BSat;
	value += Jconst;
	value +=  this->discreteGrid->muRmu0 * HnewIn;
}

template <int dimension, typename CS>
void DPC_hysteresis_model<dimension, CS>::calculate_value_in_saturation(const EigenType& HnewIn)
{
	Hnew = HnewIn;
	Hdelta = Hnew -Hold;
	value = EigenType::Zero();
	EigenType parValue;

	for (auto it = saturationSurfaces->begin(); it != saturationSurfaces->end();  ++it)
	{
		(*it).calc_value_unfrozen(Hold, Hdelta, parValue);
		value += parValue;
	}

	value *=  BSat;
	value +=  this->discreteGrid->muRmu0 * HnewIn;
	for (auto it = hysterons.begin(); it != hysterons.end(); ++it)
	{
		it->set_state_unfrozen();
	}

}

template <int dimension, typename CS>
void DPC_hysteresis_model<dimension, CS>::accept_evolution()
{
	for (int i = 0; i < numberOfHysterons; i++)
	{
		hysterons[i].accept_value(Hnew, Hdelta, Hold);
	}

	Hold = Hnew;
}

template <int dimension, typename CS>
void DPC_hysteresis_model<dimension, CS>::accept_state(EigenType HnewIn)
{
	this->Hnew = HnewIn;
	this->Hdelta = this->Hnew - this->Hold;

	for(auto it = hysterons.begin(); it != hysterons.end(); ++it)
	{
		it->update_state(Hnew, Hdelta, Hold);
	}

	this->Hold = this->Hnew;
}



template <int dimension, typename CS>
typename DPC_hysteresis_model<dimension, CS>::EigenJacob DPC_hysteresis_model<dimension, CS>::get_cont_jacobian_B_H() const
{
	EigenJacob result = EigenJacob::Zero();
	EigenJacob muMat = EigenJacob::Identity() * this->discreteGrid->muRp1mu0; // move to constants

	for (auto it = saturationSurfaces->begin(); it != saturationSurfaces->end(); ++it)
	{
		result += it->calculate_unit_vector_jacobian(Hnew);
	}
	return result * BSat + muMat;
}


template <int dimension, typename CS>
void DPC_hysteresis_model<dimension, CS>::precalculate_data(std::vector<int>& FrozenIDs, std::vector<int>& UnFrozenIDs,
														std::vector<bool> surfaceFrozen, EigenType& Jconst, EigenType HPoint)
{
	Jconst = EigenType::Zero();

	for (std::size_t i = 0; i < this->hysterons.size(); i++)
	{

		EigenType tempVal = hysterons[i].update_value(HPoint, HPoint, EigenType::Zero());
		//hysterons[i].accept_value();
		if (hysterons[i].newFrozen)
		{
			FrozenIDs.push_back(int(i));
			surfaceFrozen.push_back(true);
		}
		else
		{
			UnFrozenIDs.push_back(int(i));
			Jconst += tempVal;
			surfaceFrozen.push_back(false);
		}
		hysterons[i].reset_value(EigenType::Zero());
	}
	Jconst *=BSat;

}

template <int dimension, typename CS>
void DPC_hysteresis_model<dimension, CS>::reset_all_hysterons()
{
	this->Hnew = EigenType::Zero();
	this->Hold = EigenType::Zero();
	this->Hdelta = EigenType::Zero();
	this->value= EigenType::Zero();
	for (auto & it :hysterons )
	{
		it.reset_value(EigenType::Zero());
	}
}

template <int dimension, typename CS>
void DPC_hysteresis_model<dimension, CS>::gather_unit_magnetizations(EigenType HnewIn, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& result)
{
	this->Hnew = HnewIn;
	this->Hdelta = this->Hnew -this->Hold;
	this->value = EigenType::Zero();
	EigenType value;

	for (std::size_t i = 0; i < this->hysterons.size(); i++)
	{
		value = this->hysterons[i].update_value(this->Hnew,  this->Hdelta, this->Hold);
		result.block(0,i,2,1) = value;
	}
	this->accept_evolution();
}


template class DPC_hysteresis_model<1>;
template class DPC_hysteresis_model<2>;
template class DPC_hysteresis_model<3>;

template class DPC_hysteresis_model<1, CS_sphere<1>>;
template class DPC_hysteresis_model<2, CS_sphere<2>>;
template class DPC_hysteresis_model<3, CS_sphere<3>>;
