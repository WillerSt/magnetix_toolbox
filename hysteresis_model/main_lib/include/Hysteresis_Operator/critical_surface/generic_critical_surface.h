// SPDX-FileCopyrightText: 2025 Stephan Willerich
// SPDX-License-Identifier: MIT License

#ifndef HYSTERESIS_OPERATOR_CRITICAL_SURFACE_GENERIC_CRITICAL_SURFACE_H_
#define HYSTERESIS_OPERATOR_CRITICAL_SURFACE_GENERIC_CRITICAL_SURFACE_H_

#include <Eigen/Core>

template <int dimension>
class generic_critical_surface
{
public:
	using EigenType = Eigen::Matrix<double, dimension, 1>;
	using EigenJacob =  Eigen::Matrix<double, dimension, dimension>;

public:
	const double density;
	const unsigned int surfaceID;
	const char shapeIdentifier;
	const EigenType position;
	const double radius;

public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF(true)

public:
	generic_critical_surface(
			const double& densityIn,
			const unsigned int& surfaceIDIn,
			const char& identifierIn,
			const EigenType& positionIn,
			const double & radiusIn)
		:density(densityIn), surfaceID(surfaceIDIn), shapeIdentifier(identifierIn), position(positionIn), radius(radiusIn){}

	generic_critical_surface(
			const generic_critical_surface& surfIn)
		:density(surfIn.density), surfaceID(surfIn.surfaceID), shapeIdentifier(surfIn.shapeIdentifier), position(surfIn.position), radius(surfIn.radius){}

	inline virtual ~generic_critical_surface(){}

	bool operator == (
			const generic_critical_surface<dimension>& shapeComp) const
	{
		if (shapeIdentifier != shapeComp.shapeIdentifier)
		{
			return false;
		}
		return equality_check(shapeComp);
	}

	virtual bool form_isotropy_pair(const generic_critical_surface<dimension>& shapeComp) const
	{
		return  (( std::abs(this->position.norm() - shapeComp.position.norm()) < 1e-2) && this->same_shape_as(shapeComp));
	}

	virtual bool same_shape_as (
			const generic_critical_surface<dimension>& shapeComp) const
	{
		if (shapeIdentifier != shapeComp.shapeIdentifier)
		{
			return false;
		}
		return this->shape_check(shapeComp);
	}

	virtual bool same_position_as (
			const generic_critical_surface<dimension>& shapeComp) const
	{
		return position_check(shapeComp);
	}

	virtual bool point_symmetric_to (
			const generic_critical_surface<dimension>& shapeComp) const
	{
		return (this->symmetry_check(shapeComp)) && this->shape_check(shapeComp);
	}

	virtual bool inside(
			const EigenType& H) const = 0;

	virtual void calc_value_frozen(
			const EigenType& Hold, const EigenType& deltaH, EigenType& value) const = 0;

	virtual void calc_value_unfrozen(
			const EigenType& Hold,const EigenType& deltaH, EigenType& value) const = 0;

	virtual void calc_value_initialization(
			const EigenType& Hold, EigenType& value) const = 0;

	virtual void display_info() const = 0;

	virtual unsigned int get_surfaceID() const
	{
		return surfaceID;
	}


	virtual double get_max_radius() const = 0;

	virtual EigenJacob calculate_unit_vector_jacobian(
			const EigenType& Hnew) const = 0;

private:
	virtual bool equality_check(
			const generic_critical_surface<dimension>& surf2) const = 0;

	virtual bool shape_check(
			const generic_critical_surface<dimension>& surf2) const = 0;

	virtual bool position_check(
			const generic_critical_surface<dimension>& surf2) const = 0;

	virtual bool symmetry_check(
			const generic_critical_surface<dimension>& surf2) const = 0;

};



#endif /* HYSTERESIS_OPERATOR_CRITICAL_SURFACE_GENERIC_CRITICAL_SURFACE_H_ */
