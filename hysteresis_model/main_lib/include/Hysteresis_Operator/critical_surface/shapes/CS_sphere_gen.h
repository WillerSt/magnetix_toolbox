/*
 * sphere.h
 *
 *  Created on: 03.03.2018
 *      Author: Stephan Willerich
 */
#ifndef HYSTERESIS_OPERATOR_CRITICAL_SURFACE_SHAPES_CS_sphere_gen_GEN_H_
#define HYSTERESIS_OPERATOR_CRITICAL_SURFACE_SHAPES_CS_sphere_gen_GEN_H_

#include "Hysteresis_Operator/critical_surface/generic_critical_surface.h"
#include <cmath>


template<int dimension>
class CS_sphere_gen : public generic_critical_surface<dimension>
{

public:
	// make typedefs available
	typedef typename generic_critical_surface<dimension>::EigenType EigenType;
	typedef typename generic_critical_surface<dimension>::EigenJacob EigenJacob;

private:

	// Parameters characterizing sphere
	const double radSq;

public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF(true)

public:
	CS_sphere_gen(const EigenType & positionIn, const double& radIn, const double& densityIn, const unsigned int& surfaceIDIn)
	: generic_critical_surface<dimension>(densityIn, surfaceIDIn, 's', positionIn, radIn),
	 radSq(this->radius*this->radius){
	};

	CS_sphere_gen(const CS_sphere_gen & sphereIn)
	: generic_critical_surface<dimension>(sphereIn),radSq(this->radius*this->radius)
	  {};

	//bool operator == (const CS_sphere_gen<dimension>& sphereComp) const;
	//bool same_shape_as (const CS_sphere_gen<dimension>& shapeComp) const;
	inline bool inside(const EigenType& H) const
	{return ((H-this->position).squaredNorm()<radSq);}
	//bool same_position_as (const CS_sphere_gen<dimension>& surf2) const;

	inline void calc_value_frozen(const EigenType& Hold, const EigenType& deltaH, EigenType& value) const;

	inline void calc_value_unfrozen(const EigenType& Hold,const EigenType& deltaH, EigenType& value) const;

	void calc_value_initialization(const EigenType& Hold, EigenType& value) const;

	void display_info() const;

	EigenJacob calculate_unit_vector_jacobian(const EigenType& Hnew) const;

	double get_max_radius() const
	{
		return this->radius;
	}


private:
	inline double calc_tau(const EigenType& Hold, const EigenType& deltaH) const;
	inline void update_direction(const EigenType& Hold, const EigenType& deltaH, EigenType& value, double tau) const;
	inline bool equality_check(const generic_critical_surface<dimension>& sphereComp) const
	{
		return (
				(this->position == dynamic_cast<const CS_sphere_gen<dimension>&>(sphereComp).position)
				&& (this->radius == dynamic_cast<const CS_sphere_gen<dimension>&>(sphereComp).radius));
	}
	inline bool shape_check(const generic_critical_surface<dimension>& sphereComp) const
	{
		return (std::abs((this->radius - dynamic_cast<const CS_sphere_gen<dimension>&>(sphereComp).radius)) < 1e-10);
	}
	inline bool position_check(const generic_critical_surface<dimension>& sphereComp) const
	{
		return ((this->position - dynamic_cast<const CS_sphere_gen<dimension>&>(sphereComp).position).norm() < 1e-6);
	}
	inline bool symmetry_check(const generic_critical_surface<dimension>& sphereComp) const
	{
		return ((this->position + dynamic_cast<const CS_sphere_gen<dimension>&>(sphereComp).position).norm() < 1e-6);
	}
};
/*
template <int dimension>
bool CS_sphere_gen<dimension>::inside(const EigenType& H) const
{
	return ((H-this->position).squaredNorm()<radSq);
}
*/

template <int dimension>
void CS_sphere_gen<dimension>:: calc_value_frozen(const EigenType& Hold, const EigenType& deltaH, EigenType& value) const
{
 	double tau;
 	tau = calc_tau(Hold, deltaH);
 	update_direction(Hold, deltaH, value, tau);
}

template <int dimension>
void CS_sphere_gen<dimension>::calc_value_initialization(const EigenType& Hold, EigenType& value) const
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
			}
			else if (this->position[0] < -1e-6)
			{
				value <<   dx, -this->position[1];
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
			value = value.normalized();
		}
	}
	else
	{
		value = diff / n ;
	}

	value *= this->density;
}
template <int dimension>
double CS_sphere_gen<dimension>::calc_tau(const EigenType& Hold, const EigenType& deltaH) const
{
 	double p, root;
 	double n = deltaH.squaredNorm();
 	double nInv = 1/n;

 	if (n < 1e-12 && n > -1e-12)
 	{
 		return 1;
 	}

 	EigenType diff = Hold - this->position;

 	p = 2 * deltaH.dot(diff) * nInv;

 	root = 0.25 * p * p - (diff.squaredNorm() - this->radius * this->radius)  *  nInv;

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
void CS_sphere_gen<dimension>::calc_value_unfrozen(const EigenType& Hold,const EigenType& deltaH, EigenType& value) const
{
 	double tau = 1;
 	update_direction(Hold, deltaH, value, tau);
}

template<int dimension>
void CS_sphere_gen<dimension>::update_direction(const EigenType& Hold, const EigenType& deltaH, EigenType& value, double tau) const
{
	EigenType diff = Hold + tau * deltaH - this->position;


	double n = diff.norm();
 	// wenn n = 0 befindet sich H genau im Mittelpunkt (geht nur bei r = 0)
	if (n < 1e-12 && n > -1e-12)
 	{
 		return;
 	}

	value = diff / n * this->density;
}



#endif /* HYSTERESIS_OPERATOR_CRITICAL_SURFACE_SHAPES_CS_sphere_gen_GEN_H_ */
