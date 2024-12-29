/*
 * elementary_hysteron.h
 *
 *  Created on: 10.05.2017
 *      Author: Stephan Willerich
 */

#ifndef HEADERS_HYSTERESIS_OPERATOR_ELEMENTARY_HYSTERON_H_
#define HEADERS_HYSTERESIS_OPERATOR_ELEMENTARY_HYSTERON_H_

#include "Hysteresis_Operator/critical_surface/generic_critical_surface.h"
#include <memory>

template <int dimension, typename CS = generic_critical_surface<dimension>>
class elementary_hysteron{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF(true)

	using EigenType = typename Eigen::Matrix<double, dimension, 1>;
	using EigenJacob =  Eigen::Matrix<double, dimension, dimension>;

public:
	// variables representing the memory
	bool frozen, newFrozen;
	EigenType oldValue; //, value;

	// all behavioral information is managed within the sphere
	const std::shared_ptr<const CS> sphere;

private:
	inline void calculate_value_unfrozen();
	inline void calculate_value_frozen();

public:
	elementary_hysteron(const EigenType& initHnew, std::shared_ptr<const CS> CritSurface);

	EigenType update_value(
			const EigenType& Hnew, const EigenType& Hdelta, const EigenType& Hold);

	EigenType update_frozen_value(
			const EigenType& Hnew, const EigenType& Hdelta, const EigenType& Hold);

	void update_state(
			const EigenType& Hnew, const EigenType& Hdelta, const EigenType& Hold);

	void accept_value(
			const EigenType& Hnew, const EigenType& Hdelta, const EigenType& Hold);

	void set_state_unfrozen();

	void display_info() const;

	void reset_value(const EigenType& Hnew);

	// inline EigenType get_value() const;

	EigenJacob get_partial_jacobian_B_H(const EigenType& Hnew) const;
	// bool get_state() const;

	unsigned int get_associated_surfaceID() const;

	std::shared_ptr<const CS> get_pointer_to_associated_surface() const;



};

#include "Hysteresis_Operator/elementary_hysteron_inline.h"

#endif /* HEADERS_HYSTERESIS_OPERATOR_ELEMENTARY_HYSTERON_H_ */
