/*
 * elementary_hysteron_inline.h
 *
 *  Created on: 30.03.2018
 *      Author: Stephan Willerich
 */

#ifndef HYSTERESIS_OPERATOR_ELEMENTARY_HYSTERON_INLINE_H_
#define HYSTERESIS_OPERATOR_ELEMENTARY_HYSTERON_INLINE_H_

template <int dimension, typename CS>
inline void elementary_hysteron<dimension, CS>::set_state_unfrozen()
{
	this->newFrozen = false;
}

template <int dimension, typename CS>
inline void elementary_hysteron<dimension, CS>::update_state(const EigenType& Hnew, const EigenType& Hdelta, const EigenType& Hold)
{
	newFrozen = sphere->inside(Hnew);
	if (!frozen && newFrozen)
	{
		sphere->calc_value_frozen(Hold, Hdelta, oldValue);
	}

	frozen = newFrozen;

}

// replace old value by new value before next timestep
template <int dimension, typename CS>
inline void elementary_hysteron<dimension, CS>::accept_value(const EigenType& Hnew, const EigenType& Hdelta, const EigenType& Hold)
{
	//newFrozen = sphere->inside(Hnew);
	if (!frozen && newFrozen)
	{
		sphere->calc_value_frozen(Hold, Hdelta, oldValue);
	}

	frozen = newFrozen;
}

template<int dimension, typename CS>
inline typename elementary_hysteron<dimension, CS>::EigenJacob elementary_hysteron<dimension, CS>::get_partial_jacobian_B_H(const EigenType& Hnew) const
{
	return this->sphere->calculate_unit_vector_jacobian(Hnew);
}

template<int dimension, typename CS>
inline typename elementary_hysteron<dimension, CS>::EigenType elementary_hysteron<dimension, CS>::update_frozen_value( const EigenType& Hnew , const EigenType& Hdelta, const EigenType& Hold)
{
	newFrozen = true;

	if (frozen)
	{
		//this->value = this->oldValue;

		return oldValue;



	}


	else if (!frozen)
	{
		EigenType newVal;
		this->sphere->calc_value_frozen(Hold, Hdelta, newVal);
		//this->value = newVal;
		return newVal;
	}
	return EigenType::Zero();
}

template<int dimension, typename CS>
inline typename elementary_hysteron<dimension, CS>::EigenType elementary_hysteron<dimension, CS>::update_value( const EigenType& Hnew , const EigenType& Hdelta, const EigenType& Hold)
{
	newFrozen = sphere->inside(Hnew);

	if (frozen && newFrozen)
	{
		//this->value = this->oldValue;
		return oldValue;
	}


	else
	{
		EigenType tempVal = EigenType::Zero();
		if (!frozen && newFrozen)
		{
			this->sphere->calc_value_frozen(Hold, Hdelta, tempVal);
		}
		else
		{
			this->sphere->calc_value_unfrozen(Hold, Hdelta, tempVal);
		}

		//value = tempVal;
		return tempVal;
	}
}





#endif /* HYSTERESIS_OPERATOR_ELEMENTARY_HYSTERON_INLINE_H_ */
