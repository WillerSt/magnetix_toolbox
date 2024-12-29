/*
 * CS_sphere.h
 *
 *  Created on: 28.03.2018
 *      Author: Stephan Willerich
 */

#ifndef HYSTERESIS_OPERATOR_CRITICAL_SURFACE_SHAPES_CS_sphere_H_
#define HYSTERESIS_OPERATOR_CRITICAL_SURFACE_SHAPES_CS_sphere_H_

#include <cmath>
#include <Eigen/Core>

template<int dimension>
class CS_sphere
{

public:
	// make typedefs available
	typedef Eigen::Matrix<double, dimension, 1> EigenType;
	typedef Eigen::Matrix<double, dimension, dimension> EigenJacob;


public:
	const double density;
	const unsigned int surfaceID;
	static const char shapeIdentifier = 's';

	// Parameters characterizing sphere
	const EigenType position;
	const double radius, radSq;

public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF(true)



public:
	CS_sphere(const EigenType & positionIn, const double& radIn, const double& densityIn, const unsigned int& surfaceIDIn)
	:density(densityIn), surfaceID(surfaceIDIn),  position(positionIn), radius(radIn), radSq(radius*radius){};

	//TODO: GET RID, needed for xml- grid constructor, to be able to compile with continuous surfaces
	CS_sphere(const EigenType & positionIn, const double& radIn, const double& crap1, const double& crap2, const double& densityIn, const unsigned int& surfaceIDIn)
	:density(densityIn), surfaceID(surfaceIDIn),  position(positionIn), radius(radIn), radSq(radius*radius){};

	bool inside(const EigenType& H) const;

	bool inside(const double& distance)
	{
		return distance < radius;
	}

	bool inside_sq(const double& distance)
	{
		return distance < radSq;
	}

	double distance(const EigenType& H) const
	{
		return (H-this->position).norm();
	}

	double distanceSq(const EigenType & H) const
	{
		return (H-this->position).squaredNorm();
	}

	void distance_and_direction(const EigenType& H, EigenType& vec, double& distance) const;

	void calc_value_frozen(const EigenType& Hold, const EigenType& deltaH, EigenType& value) const;

	void calc_value_unfrozen(const EigenType& Hold,const EigenType& deltaH, EigenType& value) const;

	void calc_value_initialization(const EigenType& Hold, EigenType& value) const;

	void display_info() const;

	EigenJacob calculate_unit_vector_jacobian(const EigenType& Hnew) const;

	unsigned int get_surfaceID() const {return surfaceID;}

	bool operator == (const CS_sphere<dimension>& shapeComp) const;

	bool same_shape_as (const CS_sphere<dimension>& shapeComp) const;

	bool same_position_as (const CS_sphere<dimension>& shapeComp) const
	{return position_check(shapeComp);}

	bool point_symmetric_to (const CS_sphere& sphereComp) const
	{return (this->symmetry_check(sphereComp)) && this->shape_check(sphereComp);}

	bool form_isotropy_pair(const CS_sphere<dimension>& shapeComp) const
	{
		return  ( (this->position.norm() - shapeComp.position.norm()) < 1e-6) && same_shape_as(shapeComp);
	}

	double get_max_radius() const
	{
		return this->radius;
	}


protected:
	inline double calc_tau(const EigenType& Hold, const EigenType& deltaH) const;
	inline void update_direction(const EigenType& Hold, const EigenType& deltaH, EigenType& value, double tau) const;

	bool equality_check(const CS_sphere<dimension>& sphereComp) const;
	bool shape_check(const CS_sphere& sphereComp) const;
	bool position_check(const CS_sphere& sphereComp) const;
	bool symmetry_check(const CS_sphere& sphereComp) const
	{return ((this->position + sphereComp.position).norm() < 1e-6);}
};

#include "Hysteresis_Operator/critical_surface/shapes/CS_sphere_inline.h"


#endif /* HYSTERESIS_OPERATOR_CRITICAL_SURFACE_SHAPES_CS_sphere_H_ */
