/*
 * critical_surface_inline.h
 *
 *  Created on: 30.03.2018
 *      Author: Stephan Willerich
 */

#ifndef HYSTERESIS_OPERATOR_CRITICAL_SURFACE_INLINE_H_
#define HYSTERESIS_OPERATOR_CRITICAL_SURFACE_INLINE_H_


// Calculate value in frozen state
template <int dimension>
inline void critical_surface<dimension>::calc_value_frozen(const EigenType& Hold,
		 const EigenType& deltaH , EigenType& value) const
{
	double tau;
	tau = calc_tau(Hold, deltaH);
	update_direction(Hold, deltaH, value, tau);
}


template <int dimension>
inline void critical_surface<dimension>::calc_value_unfrozen(const EigenType& Hold,
		 const EigenType& deltaH , EigenType& value) const
{
	double tau = 1;
	update_direction(Hold, deltaH, value, tau);
}



// Calculate Parameter tau, equivalent to the cutting point between a line (evolution of H) and a circle (the sphere)
template <int dimension>
inline double critical_surface<dimension>::calc_tau(const EigenType& Hold, const EigenType& deltaH) const
{
	double p, root;
	double n = deltaH.squaredNorm();
	double nInv = 1/n;

	if (n < 1e-12 && n > -1e-12)
	{
		return 1;
	}

	EigenType diff = Hold - position;

	p = 2 * deltaH.dot(diff) * nInv;

	root = 0.25 * p * p - (diff.squaredNorm() - rad * rad)  *  nInv;

#ifdef CS_ENABLE_ERROR_CHECKING
	if (root < 0)
	{
		return -1;
	}
#endif

	if (root < 1e-10 && root > -1e-10)
	{
		 return -p * 0.5;
	}
	else
	{
		return -p * 0.5 - std::sqrt(root);
	}

}


// Update direction in accordance to tau

template <int dimension>
inline void critical_surface<dimension>::update_direction(const EigenType& Hold,
		 const EigenType& deltaH, EigenType& value, double tau) const
{
	value = Hold + tau * deltaH - position;


	double n = value.norm();
	// wenn n = 0 befindet sich H genau im Mittelpunkt (geht nur bei r = 0)
	if (n < 1e-12 && n > -1e-12)
	{
		return;
	}

	value = value  * (density / n);
}





#endif /* HYSTERESIS_OPERATOR_CRITICAL_SURFACE_INLINE_H_ */
