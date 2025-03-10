// SPDX-FileCopyrightText: 2025 Stephan Willerich
// SPDX-License-Identifier: MIT License

#include <Hysteresis_Operator/critical_surface/shapes/CS_sphere.h>
#include <iomanip>
#include <iostream>

template <int dimension>
void CS_sphere<dimension>::display_info() const
{
	std::cout << "Spherical CS No. " << this->surfaceID <<" of dimension "<< dimension << " at position [" << this->position.transpose()
			<< std::fixed << std::setprecision(8) << "]" << " with radius " << std::fixed << std::setprecision(2)
			<< this->radius <<  " and density "<< std::fixed << std::setprecision(8)
			<< this->density <<";" << std::endl;
}

template<>

typename CS_sphere<1>::EigenJacob CS_sphere<1>::calculate_unit_vector_jacobian(const EigenType& Hnew) const
 {
	return EigenJacob();
 }

template<>
typename CS_sphere<2>::EigenJacob CS_sphere<2>::calculate_unit_vector_jacobian(const EigenType& Hnew) const
 {
	 double uxx, uxy, uyx, uyy;
	 EigenJacob result;
	 EigenType diff = Hnew - this->position;
	 double diffNorm = diff.norm();
	 if (diffNorm > 1e-6)
	 {
		 uxx = (diffNorm - std::pow(diff[0],2) / diffNorm) / std::pow(diffNorm,2);
		 uxy = (-diff[0] * diff[1] / diffNorm) / std::pow(diffNorm,2);
		 uyx =uxy;
		 uyy = (diffNorm - std::pow(diff[1],2) / diffNorm) / std::pow(diffNorm,2);
	 }
	 else
	 {
		 uxx = 0;
		 uxy = 0;
		 uyx = 0;
		 uyy = 0;
	 }

	 result << uxx, uxy, uyx, uyy;

	 result *= this->density;
	 return result;
 }

template<>
typename CS_sphere<3>::EigenJacob CS_sphere<3>::calculate_unit_vector_jacobian(const EigenType& Hnew) const
 {
	return EigenJacob();
 }

template class CS_sphere<1>;
template class CS_sphere<2>;
template class CS_sphere<3>;
