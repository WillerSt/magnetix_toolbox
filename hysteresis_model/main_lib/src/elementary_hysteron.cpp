/*
 * elementary_hysteron.cpp
 *
 *  Created on: 10.05.2017
 *      Author: Stephan Willerich
 */
#include "Hysteresis_Operator/elementary_hysteron.h"
#include "Hysteresis_Operator/critical_surface/known_shapes.h"
#include "iostream"


// Constructor sets all variables
template <int dimension, typename CS>
elementary_hysteron<dimension, CS>::elementary_hysteron(const EigenType& Hnew,
		std::shared_ptr<const CS> CritSurface)
: frozen(false), newFrozen(false), sphere(CritSurface)
{
	reset_value(Hnew);
}


// display information of the sphere
template<int dimension, typename CS>
void elementary_hysteron<dimension, CS>::display_info() const
{
	sphere->display_info();
}

template<int dimension, typename CS>
void elementary_hysteron<dimension, CS>::reset_value(const EigenType& Hnew)
{
	frozen = (sphere)->inside(Hnew);
	newFrozen = frozen;

	sphere->calc_value_initialization(Hnew, this->oldValue);

	//this->oldValue = this->value;

}



template <int dimension, typename CS>
unsigned int elementary_hysteron<dimension, CS>::get_associated_surfaceID() const
{
	return sphere->get_surfaceID();
}

template <int dimension, typename CS>
std::shared_ptr<const CS> elementary_hysteron<dimension, CS>::get_pointer_to_associated_surface() const
{
	return this->sphere;
}

template class elementary_hysteron<1>;
template class elementary_hysteron<2>;
template class elementary_hysteron<3>;

template class elementary_hysteron<1, CS_sphere<1>>;
template class elementary_hysteron<2, CS_sphere<2>>;
template class elementary_hysteron<3, CS_sphere<3>>;
