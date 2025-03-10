// SPDX-FileCopyrightText: 2025 Stephan Willerich
// SPDX-License-Identifier: MIT License

#include <Hysteresis_Operator/malloc_grid.h>
#include "Hysteresis_Operator/critical_surface/known_shapes.h"
#include  "hystLib_constants.h"
#include <iostream>

template <int dimension, typename CS>
malloc_grid<dimension, CS>::malloc_grid(dpc_grid_constructor<dimension, CS>& GridConstructor, double BsatIn)
:pointCollection(GridConstructor.memoryAllocationGrid),
 surfaces(GridConstructor.surfaces),
 saturationSurfaces(construct_saturation_surfaces()),
 Bsat(BsatIn), numberOfSurfaces(surfaces->size()),
 muR(GridConstructor.get_lin_part()*this->Bsat),
 muRmu0(hystLib::constants::mu0 * this->muR),
 muRp1mu0(hystLib::constants::mu0 * (this->muR+1))
{
	bool orderError = false;
	for (std::size_t i = 0; i < saturationSurfaces->size(); i++)
	{
		if ( ((*saturationSurfaces)[i].position - pointCollection[i].coord).norm() > 1e-10)
		{
			orderError = true;
			break;
		}
	}

	if (orderError)
	{
		std::cout << "WARNING: Order of point collection and unfrozen surfaces in malloc grid does not match\n";
	}
	else
	{
		// std::cout << "INFO: Ordering of vectors in malloc grid is consistent \n";
	}

	orderError = false;
	bool breakOuter = false;

	for (std::size_t i = 0; i < pointCollection.size(); i++)
	{
		for(std::size_t j = pointCollection[i].surfIndBegin; j < pointCollection[i].surfIndEnd; j++ )
		{
			if ((*surfaces)[j]->radius > (*surfaces)[j+1]->radius)
			{
				std::cout << j << " " << j+1;
				breakOuter = true;
				orderError = true;
				break;
			}
		}
		if (breakOuter)
		{
			break;
		}
	}

	if (orderError)
	{
		std::cout << "WARNING: Order of critical surfaces in malloc grid does not match\n";
	}
	else
	{
		//std::cout << "INFO: Ordering of critical surfaces in malloc grid is consistent \n";
	}

	bool partialSumError = false;

	for (std::size_t i = 0; i < pointCollection.size(); i++)
	{
		if ( std::abs(pointCollection[i].partialDensSums[pointCollection[i].partialDensSums.size()-1] - (*saturationSurfaces)[i].density) > 1e-14)
		{
			partialSumError = true;
			break;
		}
	}

	if (partialSumError)
	{
		std::cout << "WARNING: Partial density sums in malloc grid does not match\n";
	}
	else
	{
		// std::cout << "INFO: Partial density sums in malloc grid are consistent \n";
	}



}

template <int dimension, typename CS>
std::shared_ptr<std::vector<unfrozen_surface<dimension>>> malloc_grid<dimension, CS>::construct_saturation_surfaces()
{
	auto saturationSurfaces = std::make_shared<std::vector<unfrozen_surface<dimension>>> ();
	for(unsigned int i  = 0; i < this->pointCollection.size(); i++)
	{
		EigenType coord = pointCollection[i].coord;
		double densSum  = 0;
		double radMax = 0;
		for (unsigned int j = pointCollection[i].surfIndBegin; j <= pointCollection[i].surfIndEnd; j++)
		{
			if ((*surfaces)[j]->get_max_radius() > radMax)
			{
				radMax = (*surfaces)[j]->get_max_radius();
			}
			densSum += (*surfaces)[j]->density;
		}
		saturationSurfaces->push_back(unfrozen_surface<dimension>(coord, densSum, i, radMax)); // @suppress("Method cannot be resolved")
	}
	return saturationSurfaces;
}

template <int dimension, typename CS>
void malloc_grid<dimension, CS>::display_all_surface_infos() const
{
	for (unsigned int i = 0; i < malloc_grid::pointCollection.size(); i++)
	{
		std::cout << "Memory Grid Point No. " << i+1 << " with critical surfaces:" << std::endl;
		//malloc_grid:: pointCollection[i].display_info(2);
	}
}

template <int dimension, typename CS>
void malloc_grid<dimension, CS>::display_reduced_surface_infos() const
{
	for (unsigned int i = 0; i < malloc_grid::pointCollection.size(); i++)
	{
		std::cout << "Memory Grid Point No. " << i+1 << " with critical surfaces:" << std::endl;
		//malloc_grid:: pointCollection[i].display_info(1);
	}
}


template <int dimension, typename CS>
void malloc_grid<dimension, CS>::display_point_info(int i) const
{
	if (i< 0 || (i > int(malloc_grid::pointCollection.size()) - 1))
	{
		std::cout << "Error: Queried Index "  << i << " out of range "
				<< 0<<" - " <<int(malloc_grid::pointCollection.size()) - 1 << "."<< std::endl;
	}
	else
	{
		//malloc_grid:: pointCollection[i].display_info(2);
	}

}


template <int dimension, typename CS>
typename malloc_grid<dimension, CS>::point_information malloc_grid<dimension, CS>::access_grid_point(int i) const
{
	if (i< 0 || (i > int(malloc_grid::pointCollection.size()) - 1))
	{
		std::cout << "Error: Grid Point "  << i << " does not exist, range is: "
				<< 0<<" - " <<int(malloc_grid::pointCollection.size()) - 1 << "."<< std::endl;
	return malloc_grid::pointCollection[0];
	}
	else
	{
	return malloc_grid::pointCollection[i];
	}
}


template <int dimension, typename CS>
unsigned int malloc_grid<dimension, CS>::get_number_of_points() const
{
	return (malloc_grid::pointCollection).size();
}


template <int dimension, typename CS>
std::shared_ptr<const CS> malloc_grid<dimension, CS>::get_pointer_to_surface(unsigned int i) const
{
	return (*this->surfaces)[i];
}

template <int dimension, typename CS>
std::shared_ptr<std::vector<unfrozen_surface<dimension>>> malloc_grid<dimension, CS>::get_pointers_to_saturation_surfaces() const
{
	return saturationSurfaces;
}


template class malloc_grid<1, generic_critical_surface<1>>;
template class malloc_grid<2, generic_critical_surface<2>>;
template class malloc_grid<3, generic_critical_surface<3>>;

template class malloc_grid<1, CS_sphere<1>>;
template class malloc_grid<2, CS_sphere<2>>;
template class malloc_grid<3, CS_sphere<3>>;

template class malloc_grid<1, CS_sphere_cont<1>>;
template class malloc_grid<2, CS_sphere_cont<2>>;
template class malloc_grid<3, CS_sphere_cont<3>>;
