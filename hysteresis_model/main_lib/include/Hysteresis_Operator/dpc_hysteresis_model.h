// SPDX-FileCopyrightText: 2025 Stephan Willerich
// SPDX-License-Identifier: MIT License

#ifndef HEADERS_HYSTERESIS_OPERATOR_DPC_HYSTERESIS_MODEL_
#define HEADERS_HYSTERESIS_OPERATOR_DPC_HYSTERESIS_MODEL_

#include "Hysteresis_Operator/malloc_grid.h"
#include "Hysteresis_Operator/elementary_hysteron.h"
#include "hystLib_constants.h"
#include <memory>


template <int dimension, typename CS = generic_critical_surface<dimension>>
class DPC_hysteresis_model
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF(true)
	using EigenType = Eigen::Matrix<double, dimension, 1>;
	using EigenJacob = Eigen::Matrix<double, dimension, dimension>;
	using sphereType = CS;

protected:
	EigenType Hnew = EigenType::Zero();
	EigenType Hold = EigenType::Zero();
	EigenType Hdelta = EigenType::Zero();
	EigenType value= EigenType::Zero();

	const std::shared_ptr<const malloc_grid<dimension, CS>> discreteGrid;

	int numberOfHysterons = 0;

public:
	const double BSat;

protected:
	std::vector<elementary_hysteron<dimension, CS>> hysterons;
	std::shared_ptr<std::vector<unfrozen_surface<dimension>>> saturationSurfaces;




public:


	DPC_hysteresis_model(const std::shared_ptr<const malloc_grid<dimension, CS>>& malloc_gridIN);

	void calculate_evolution_H(const EigenType& HnewIn){calculate_evolution(HnewIn);}

	void calculate_evolution(const EigenType& HnewIn);

	void calculate_evolution(const EigenType& HnewIn, const std::vector<int>& FrozenIDs, const EigenType& Jconst);

	void calculate_value_in_saturation(const EigenType& HnewIn);

	void calculate_nonlinear_and_constant_part(const EigenType& HnewIn, EigenType& JNonLinear,
		const std::vector<int>& FrozenSurfaceIDs, const std::vector<int>& UnFrozenSurfaceIDs );

	void accept_state(EigenType HnewIn);

	void precalculate_data(	std::vector<int>& FrozenIDs, std::vector<int>& UnFrozenIDs,
							std::vector<bool> surfaceFrozen, EigenType& Jconst, EigenType HPoint);

	void accept_evolution();

	void reset_all_hysterons();

	EigenType get_H() const;

	EigenType get_M() const;

	EigenType get_B() const;

	EigenType get_J() const;

	EigenJacob get_cont_jacobian_B_H() const;

	int get_number_of_hysterons();

	void gather_unit_magnetizations(EigenType HnewIn, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& result);

protected:
	std::vector<elementary_hysteron<dimension, CS>> initialize_hysterons();
};
#include "Hysteresis_Operator/dpc_hysteresis_model_inline.h"

#endif /* HEADERS_HYSTERESIS_OPERATOR_DPC_HYSTERESIS_MODEL_ */
