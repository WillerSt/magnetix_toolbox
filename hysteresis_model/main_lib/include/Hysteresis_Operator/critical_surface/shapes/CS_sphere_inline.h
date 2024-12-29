/*
 * CS_sphere_inline.h
 *
 *  Created on: 30.03.2018
 *      Author: Stephan Willerich
 */

#ifndef HYSTERESIS_OPERATOR_CRITICAL_SURFACE_SHAPES_CS_SPHERE_INLINE_H_
#define HYSTERESIS_OPERATOR_CRITICAL_SURFACE_SHAPES_CS_SPHERE_INLINE_H_
#include <iostream>
/*
template <int dimension>
bool CS_sphere<dimension>::inside(const EigenType& H) const
{
	return ((H-this->position).squaredNorm()<radSq);
}
*/

template <int dimension>
inline bool CS_sphere<dimension>::inside(
		const EigenType& H) const
{
	return ((H-this->position).squaredNorm() < radSq);
}



template <int dimension>
inline void CS_sphere<dimension>:: calc_value_frozen(
		const EigenType& Hold, const EigenType& deltaH, EigenType& value) const
{
 	update_direction(Hold, deltaH, value, calc_tau(Hold, deltaH));
}



template <int dimension>
void CS_sphere<dimension>::calc_value_initialization(
		const EigenType& Hold, EigenType& value) const
{
	EigenType diff = Hold - this->position;
	double n = diff.norm();
	// wenn n = 0 befindet sich H genau im Mittelpunkt (geht nur bei r = 0)
	if (n < 1e-12 && n > -1e-12)
	{
		/* Gives best results with identification
		if (this->surfaceID%2 == 0)
			value << 1* this->density,0;
		else
			value << -1* this->density,0;
		*/
		value = EigenType::Zero();
		return;
	}

	if (this->inside(Hold))
	{
		EigenType Hdelta;
		EigenType dirVec;

		if (std::abs(this->position[1]) < 1e-6)
		{
			if (this->position[0] > 1e-6)
			{
				value << -1,0;


			}
			else
			{
				value << 1,0;
			}
		}
		else
		{
			double dx = sqrt(this->radius*this->radius - this->position[1]*this->position[1]);
			if (this->position[0] > 1e-6)
			{
				value << - dx, -this->position[1];
				value = value.normalized();
			}
			else if (this->position[0] < -1e-6)
			{
				value <<   dx, -this->position[1];
				value = value.normalized();
			}
			else
			{
				/*
				if (this->surfaceID%2 == 0)
					value <<   dx, -this->position[1];
				else
					value << - dx, -this->position[1];
				 */
				value = EigenType::Zero();
			}

		}
	}
	else
	{
		value = diff / n ;
	}

	value *= this->density;
}


template <int dimension>
inline void  CS_sphere<dimension>::distance_and_direction(const EigenType& H, EigenType& vec, double& distance) const
{
	vec = H - this->position;
	distance = vec.norm();
	if (distance < 1e-6)
	{
		vec = EigenType::Zero();

	}
	else
	{
		vec = vec/ distance;
	}
}


template <int dimension>
inline double CS_sphere<dimension>::calc_tau(
		const EigenType& Hold, const EigenType& deltaH) const
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

 	root = 0.25 * p * p - (diff.squaredNorm() - radius * radius)  *  nInv;

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



template<int dimension>
inline void CS_sphere<dimension>::calc_value_unfrozen(
		const EigenType& Hold,const EigenType& deltaH, EigenType& value) const
{
 	update_direction(Hold, deltaH, value, 1);
}



template<int dimension>
inline void CS_sphere<dimension>::update_direction(
		const EigenType& Hold, const EigenType& deltaH, EigenType& value, double tau) const
{
	EigenType diff = Hold + tau * deltaH - this->position;


	double n = diff.norm();
 	// wenn n = 0 befindet sich H genau im Mittelpunkt (geht nur bei r = 0)
	if (n < 1e-12 && n > -1e-12)
 	{
		value = EigenType::Zero();
 		return;
 	}

	value = diff / n * this->density;


}



template<int dimension>
inline bool CS_sphere<dimension>::operator == (
		const CS_sphere<dimension>& shapeComp) const
{
	if (shapeIdentifier != shapeComp.shapeIdentifier)
	{
		return false;
	}
	return equality_check(shapeComp);
}



template<int dimension>
inline bool CS_sphere<dimension>::same_shape_as (
		const CS_sphere<dimension>& shapeComp) const
{
	if (shapeIdentifier != shapeComp.shapeIdentifier)
	{
		return false;
	}
	return shape_check(shapeComp);
}



template<int dimension>
inline bool CS_sphere<dimension>::equality_check(
		const CS_sphere<dimension>& sphereComp) const
{
	return (
			(this->position == sphereComp.position)
			&& (this->radius == sphereComp.radius));
}



template<int dimension>
inline bool CS_sphere<dimension>::shape_check(
		const CS_sphere& sphereComp) const
{
	return (std::abs((this->radius - sphereComp.radius)) < 1e-10);
}



template<int dimension>
inline bool CS_sphere<dimension>::position_check(
		const CS_sphere& sphereComp) const
{
	return ((this->position-(sphereComp).position).norm() < 1e-6);
}



#endif /* HYSTERESIS_OPERATOR_CRITICAL_SURFACE_SHAPES_CS_SPHERE_INLINE_H_ */
