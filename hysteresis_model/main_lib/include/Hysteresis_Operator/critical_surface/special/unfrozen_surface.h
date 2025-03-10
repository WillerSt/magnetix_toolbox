// SPDX-FileCopyrightText: 2025 Stephan Willerich
// SPDX-License-Identifier: MIT License

#ifndef HYSTERESIS_OPERATOR_CRITICAL_SURFACE_SPECIAL_UNFROZEN_SURFACE_H_
#define HYSTERESIS_OPERATOR_CRITICAL_SURFACE_SPECIAL_UNFROZEN_SURFACE_H_

#include "Hysteresis_Operator/critical_surface/shapes/CS_sphere.h"

template <int dimension>
class unfrozen_surface : public CS_sphere<dimension>
{
public:
	using EigenType =  Eigen::Matrix<double, dimension, 1>;
	using EigenJacob = Eigen::Matrix<double, dimension, dimension>;

	const double radMax;
	const double radMaxSq = radMax*radMax;

	unfrozen_surface(const EigenType & positionIn, const double& densityIn, const unsigned int& surfaceIDIn, const double& radMaxIn)
	:CS_sphere<dimension>(positionIn, 0, densityIn, surfaceIDIn), radMax(radMaxIn){};
};



#endif /* HYSTERESIS_OPERATOR_CRITICAL_SURFACE_SPECIAL_UNFROZEN_SURFACE_H_ */
